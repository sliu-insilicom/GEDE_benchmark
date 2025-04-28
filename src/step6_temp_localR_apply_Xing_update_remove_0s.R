# The galaxy3 server is down, so I put a local R + VS code (doesn't fit Rstudio well)
# so use this R file for analysis - will merge code back to Rmd when server is back up


# Load the required libraries
library(AnVIL)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library("recount3")
library("limma")
library("SummarizedExperiment")  # assays function come from here
library(pheatmap)
library(DESeq2)
library(Matrix)
library(GEDE)
library(edgeR)
library(xgboost)
library(randomForest)
library(Hmisc)
library(VennDiagram)
library(eulerr)
library(patchwork)
library(biomaRt)

# src functions
setwd("/home/sliu/github/GEDE_benchmark/src")
source("GEDE_data_prep.R")
source("GEDE_simulation.R")
source("GEDE_real_analysis.R")

# change output dir
out_dir <- "/home/sliu/github/GEDE_benchmark/output/real_data_analysis"
setwd(out_dir)


####### below are scripts to run the real data analysis
set.seed(3369)

# Load the data
objs <- load("input_Recount3_BRCA_race.Rdata")
input_data <- Data_BRCA
race = "tcga.gdc_cases.demographic.race"
deg_condition = "tcga.gdc_cases.demographic.race"
cutoff_pvalue = 0.05
cutoff_log2fc = 0.5


# remove outlier:counts## Winsorization (replacing outliers by data-driven thresholds) based
## on Hampel filter
Winsor <- function(Y, nMAD=3) {
  if (is.null(nMAD)) { #no filter
    out.upper.ids <- out.lower.ids <- integer(0)
  } else {
    Y <- as.matrix(Y)
    Ymed <- rowMedians(Y, na.rm=TRUE)
    Ymad <- rowMads(Y, na.rm=TRUE)
    U <- Ymed+nMAD*Ymad; L <- Ymed-nMAD*Ymad
    out.upper.ids <- which(Y-U>0, arr.ind=TRUE)
    out.lower.ids <- which(Y-L<0, arr.ind=TRUE)
    ## replace outliers by low/upper thresholds
    for (i in 1:nrow(out.upper.ids)) {
      i2 <- out.upper.ids[i,1]; j <- out.upper.ids[i,2]
      Y[i2,j] <- U[i2]
    }
    ##
    for (i in 1:nrow(out.lower.ids)) {
      i2 <- out.lower.ids[i,1]; j <- out.lower.ids[i,2]
      Y[i2,j] <- L[i2]
    }
  }
  return(list(Y=Y, out.upper.ids=out.upper.ids, out.lower.ids=out.lower.ids, U=U, L=L))
}



######## here is the new adjustment 1/2: we remove samples with low depth
# prep data for GEDE
df_count0 = input_data$counts; m0 <- nrow(df_count0); n0 <- ncol(df_count0)
df_meta = input_data$group
temp_metadata <- df_meta[colnames(df_count0), , drop = TRUE]  # Reorder to match count data

## # this is group vector
## temp_group <- factor(temp_metadata[[deg_condition]])
temp_group <- factor(temp_metadata[[race]])

## limma does not like spaces
levels(temp_group) <- c("black", "white")

# numeric for GEDE input: white: 1; black: 0
numeric_group <- as.numeric(temp_group) - 1

## identify samples with extremely low total number of reads (in millions)
TNR <- colSums(df_count0)/1e6
o <- Winsor(t(TNR))
samples.low <- o$out.lower.ids[, 2]
## 
df_count1 <- df_count0[, -samples.low]; n <- ncol(df_count1)
temp_metadata1 <- temp_metadata[-samples.low,]
temp_group1 <- temp_group[-samples.low]
numeric_group1 <- numeric_group[-samples.low]

## x11 doesn't work for my local R
pdf("update_fig_1_histogram_of_depth.pdf", width = 8, height = 6)
hist(TNR, 51, xlab = "Total Reads", main = "")
abline(v = o$L, lty = 2)
dev.off()


######## here is the new adjustment 2/2: we remove genes with more than half 0s
## Xing 04/12/2025. Filtering genes with excessive zeros
n.zeros <- rowsums(df_count1==0)
dim(df_count1)  # 22321

df_count <- df_count1[n.zeros <= n/2,]; m <- nrow(df_count)
dim(df_count) # 17978
# 4343 genes are removed

# log2 transformation: GEDE need log data because it assumes multi-normal dist
log_df = log2(df_count+1)

## log-transformation improves normality. 
W0 <- apply(df_count, 1, function(x) as.numeric(shapiro.test(x)$statistic))
W1 <- apply(log_df, 1, function(x) as.numeric(shapiro.test(x)$statistic))

kable(rbind(Orig=summary(W0), Log=summary(W1)), caption="Summary of the W statistics of the original counts and log-transformed data. W statistics reflect the correlation between the empirical distribution and theoretical normal distribution. A large W statistic implies that the data are more normally distributed.")

## randomly select a few genes for visualization
set.seed(123)
idx <- sample(m, 6)
Y0 <- df_count[idx,]; Y1 <- log_df[idx, ]
# Open the PDF device
pdf("update_fig_2_visualize_few_genes_normality.pdf", width = 8, height = 12)  # Tall page because 6x2 layout
par(mfrow = c(6, 2))
for (i in 1:length(idx)) {
  g <- idx[i]
  hist(Y0[i, ], breaks = 21, xlab = "Counts",
       main = paste0("Gene ", g, ", W=", round(W0[g], 3)))
  hist(Y1[i, ], breaks = 21, xlab = "Log-expression",
       main = paste0("Gene ", g, ", W=", round(W1[g], 3)))
}
dev.off()


## run GEDE on filtered data
t1 <- system.time(gede_output <- GEDE(t(log_df), X=numeric_group1, K.method="vprop", HD=TRUE))

# GEDE improved expressions
new_df = transpose(gede_output$Ystar)
rownames(new_df) <- rownames(log_df); colnames(new_df) <- colnames(new_df)

## identify outliers based on marginal distribution of log_df, using a
## liberal 2 MAD cutoff
outliers2 <- as.data.frame(Hampel(t(log_df), nMAD=2, arr.ind=TRUE))
colnames(outliers2) <- c("sample", "gene")
## 6.3% of data are marked as potential outliers
round(nrow(outliers2)/length(log_df)*100, 1)

# transfer back to count by 2**x-1, floored at 0
new_df_count = round(pmax(2^new_df-1, 0))

# Open the PDF device
pdf("update_fig_3_visualize_genes_before_after_GEDE.pdf", width = 8, height = 10)  # adjust width/height if you want
# Set the plotting layout
par(mfrow = c(3, 2))  # 3 rows, 2 columns of plots per page
# Your plotting loop
for (i in idx) {
  yi <- new_df[i, ]
  xi <- log_df[i, ]
  xy <- range(c(xi, yi))
  plot(yi ~ xi, xlab = "Before GEDE", ylab = "After GEDE", xlim = xy, ylim = xy, main = i)
  out <- subset(outliers2, gene == i)[["sample"]]
  points(yi[out] ~ xi[out], pch = 19, col = 2)
  abline(0, 1)
}
dev.off()


### LFC-var-comparison
# LFC (log-fold-changes) before/after GEDE, using log-data
n0 <- sum(numeric_group1==0) #191 black subjects
n1 <- sum(numeric_group1)    #298 white subjects
x <- ifelse(numeric_group1==0, -1/n0, 1/n1)
LFC.orig <- drop(log_df%*%x); LFC.gede <- drop(new_df%*%x)

## Also compare the pooled variance
X <- cbind(1, numeric_group1); Hat <- X%*%solve(crossprod(X))%*%t(X)
yhat0 <- log_df%*%Hat; yhat1 <- new_df%*%Hat
res0 <- log_df-yhat0; res1 <- new_df-yhat1
var0 <- rowSums(res0^2)/(n-2); var1 <- rowSums(res1^2)/(n-2)

# Open the PDF device
pdf("update_fig_4_GEDE_LFC_variance_plots.pdf", width = 10, height = 5)  # Wider page for 1x2 layout
# Set up the plotting layout: 1 row, 2 columns
par(mfrow = c(1, 2))
# First plot: LFC before vs after GEDE
xy <- range(c(LFC.orig, LFC.gede))
plot(LFC.gede ~ LFC.orig, xlim = xy, ylim = xy, main = "LFC")
abline(0, 1)
# Second plot: Variance before vs after GEDE
xy <- range(c(var0, var1))
plot(var1 ~ var0, xlim = xy, ylim = xy,
     xlab = "Before GEDE", ylab = "After GEDE",
     main = "Pooled sample variance")
abline(0, 1)
dev.off()





###################### Now we can come back to my own analysis
# I didn't use Xing's ORACLE output in L244
# df_count is df after gene removel
dim(df_count)

# prep data for GEDE
df_meta = input_data$group
temp_metadata <- df_meta[colnames(df_count), , drop = TRUE]  # Reorder to match count data
# this is group vector
temp_group <- factor(temp_metadata[[deg_condition]])
levels(temp_group) <- gsub(" ", "_", levels(temp_group))  # replace space by _ for the limma contrast
# numeric for GEDE input
numeric_group <- as.numeric(temp_group) - 1

# GEDE needs transposed input: gene * sample
transpose_df_count <- transpose(df_count)
# Restore row and column names
rownames(transpose_df_count) <- colnames(df_count)  # Old column names → new row names
colnames(transpose_df_count) <- rownames(df_count)  # Old row names → new column names


# log2 transformation: GEDE need log data because it assumes multi-normal dist
log_df = log2(transpose_df_count+1)
normality_results <- apply(log_df, 2, function(x) {
  if (length(unique(x)) > 1) {  # Ensure column is not constant
    return(shapiro.test(x)$p.value)
  } 
})
normality_results <- unlist(normality_results)
dim(log_df)

# run gede based on log df and then transfer it back
gede_output <- GEDE(log_df, X=numeric_group, K.method="vprop", HD=TRUE)

# check output (transpose back)
new_df = transpose(gede_output$Ystar)
rownames(new_df) <- colnames(gede_output$Ystar)  # Old column names → new row names
colnames(new_df) <- rownames(gede_output$Ystar)  # Old row names → new column names

# transfer back to count by 2**x
new_count_df = round(pmax(2^new_df-1, 0))

# change of depth
# change of depth (colsums)
dim(df_count)
dim(new_count_df)
old_depth = colsums(df_count)
new_depth = colsums(new_count_df)
new_ratio = new_depth / old_depth

# it's not suprise to see the overall depth is decrease - we are removing outliers, combine with the LFC change by group
pdf("update_fig_5_GEDE_depth_ratio_histogram.pdf", width = 7, height = 5)  # Reasonable size for a single histogram
# Your histogram code
hist(new_ratio, 
     breaks = 100,         # Number of bins
     col = "blue",         # Bar fill color
     border = "black",     # Border color
     main = "Histogram of depth ratios", 
     xlab = "Ratio", 
     xlim = c(0, 2),       # Limit x-axis from 0 to 2
     ylab = "Frequency")
dev.off()


# fold changes before/after GEDE
old_means <- sapply(unique(numeric_group), function(g) {
  rowMeans(df_count[, numeric_group == g, drop = FALSE])  # Compute mean for each row within the group
})
old_fc <- old_means[,1]  / (old_means[,2]+1)

new_means <- sapply(unique(numeric_group), function(g) {
  rowMeans(new_count_df[, numeric_group == g, drop = FALSE]) 
})
new_fc <- new_means[,1] / (new_means[,2]+1)


# plot differences
plot_old_fc = old_fc
plot_new_fc = new_fc
# make plot
plot_old_fc[plot_old_fc > 10] <- 10
plot_old_fc[plot_old_fc < -10] <- -10
plot_new_fc[plot_new_fc > 10] <- 10
plot_new_fc[plot_new_fc < -10] <- -10

# Combine into a dataframe for boxplot
df <- data.frame(
  FoldChange = c(plot_old_fc, plot_new_fc),
  Group = rep(c("Old FC", "New FC"), each = length(plot_old_fc))
)
medians <- tapply(df$FoldChange, df$Group, median)

# Create boxplot
pdf("update_fig_6_GEDE_foldchange_boxplot.pdf", width = 7, height = 6)  # Good size for a boxplot

boxplot(FoldChange ~ Group, data = df, 
        main = "Comparison of Old and New Fold Change (Clipped at ±10)", 
        ylab = "Fold Change",
        col = c("lightblue", "lightgreen"),  # Box colors
        border = "black",                    # Outline color
        notch = TRUE)                         # Notched boxplot for median confidence

dev.off()

# trim the metadata for matching rows and cols
df_meta <- df_meta[rownames(df_meta) %in% colnames(df_count), ]
dim(df_meta)


# here we winsorize the df_meta as well 
log_df = log2(df_count+1)
dim(log_df)
o <- Winsor(log_df, nMAD=3)
## around 1.3% of data were winsorized
(nrow(o$out.upper.ids)+nrow(o$out.lower.ids))/length(log_df)

Yw <- o$Y #this is the set of winsorized log-expressions
## the winsorized count data
df_count.w <- round(pmax(2^Yw-1, 0))
dim(df_count.w)

# use counts = df_count for un-winsorized data
new_input = list(group = df_meta, counts = df_count.w, gede_counts = new_count_df)
Data_BRCA <- new_input
all(sum(rownames(Data_BRCA$group)  == colnames(Data_BRCA$gede_counts)))
save(Data_BRCA, file="updated_input_Recount3_BRCA_race.Rdata")
input_data = new_input



### Can pick up here for future analysis
# load("updated_input_Recount3_BRCA_race.Rdata")
# input_data <- Data_BRCA





### 2.1_call_DEG_for_ground_truth_on_full_data
# clear old results
rm(deg_deseq2, deg_limma, deg_edger)

deg_deseq2 <- run_DESeq2_for_df_and_meta(df_count=input_data$counts, df_meta=input_data$group, single_condition=deg_condition)

# limma needs a voom transformed data
deg_limma <- run_limma_for_df_and_meta(df_count=input_data$counts, df_meta=input_data$group, single_condition=deg_condition)

# edgeR
deg_edger <- run_edgeR_for_df_and_meta(df_count=input_data$counts, df_meta=input_data$group, single_condition=deg_condition)


# cutoff 0.05 and 0.5
filtered_deseq2 <- filter_deg(deg_deseq2, log2fc=cutoff_log2fc)
print(dim(filtered_deseq2))
filtered_limma <- filter_deg(deg_limma, log2fc = cutoff_log2fc)
print(dim(filtered_limma))
filtered_edger <- filter_deg(deg_edger, log2fc = cutoff_log2fc)
print(dim(filtered_edger))


# re-define the venn-plot function to save to PDF
compare_deg_non_proportional_figure <- function(v1, v2, v3, 
                                                 names = c("DESeq2", "limma", "EdgeR"),
                                                 pdf_file = "compare_deg_venn.pdf") {
  # Open PDF device
  pdf(pdf_file, width = 7, height = 7)  # You can adjust the size
  
  # Compute overlaps
  n12 <- length(intersect(v1, v2))  # v1 ∩ v2
  n13 <- length(intersect(v1, v3))  # v1 ∩ v3
  n23 <- length(intersect(v2, v3))  # v2 ∩ v3
  n123 <- length(Reduce(intersect, list(v1, v2, v3)))  # v1 ∩ v2 ∩ v3
  
  # Create new page
  grid.newpage()
  
  # Draw the Venn diagram
  venn.plot <- draw.triple.venn(
    area1 = length(v1),
    area2 = length(v2),
    area3 = length(v3),
    n12 = n12, 
    n13 = n13, 
    n23 = n23, 
    n123 = n123, 
    category = names,  # Method names
    fill = c("red", "blue", "green"),  # Colors
    alpha = 0.5,
    lwd = 2,
    cex = 1.5,
    cat.cex = 1.5,
    cat.col = c("red", "blue", "green")
  )
  
  grid.draw(venn.plot)
  
  # Close PDF device
  dev.off()
}


# regular Venn
compare_deg_non_proportional_figure(rownames(filtered_deseq2), rownames(filtered_limma), rownames(filtered_edger), pdf_file = "update_venn_before_GEDE.pdf")

# try the log2fc 1 plot
compare_deg_non_proportional_figure(rownames(filter_deg(deg_deseq2, log2fc=1)), rownames(filter_deg(deg_limma, log2fc = 1)), rownames(filter_deg(deg_edger, log2fc = 1)), pdf_file = "update_venn_before_GEDE_log2_cutoff1.pdf")


# based on RNA-seq data results, limm results is much smaller than the other 2, so I will use DESeq2+EdgeR overlap/union but ignore limma
out6 = list()
out6$deg_all_records_deseq2 <- deg_deseq2
out6$deg_all_records_limma <- deg_limma
out6$deg_all_records_edger <- deg_edger

out6$deg_sig_deseq2 <- filtered_deseq2
out6$deg_sig_limma <- filtered_limma
out6$deg_sig_edger <- filtered_edger


# here we still use overlaps of 2 for convenience

out6$union_deg <- union(rownames(filtered_deseq2), rownames(filtered_edger))  # this is the finding before data correction
out6$overlap_deg <- intersect(rownames(filtered_deseq2), rownames(filtered_edger)) # this serves as ground truth finding

save(out6, file="output_real_data_analysis.Rdata")



### 2.2, subsample saturation for DEG before GEDE
df_meta <- input_data$group
group_values = levels(as.factor(df_meta[[deg_condition]]))
# split group and do random sample
sample_grp_1 = rownames(df_meta[df_meta[deg_condition] == group_values[1], ])
sample_grp_2 = rownames(df_meta[df_meta[deg_condition] == group_values[2], ])

# loop DEGs by different sample rate, note that for the 200 vs 300 RNA-seq case, 20% can recover 87.5% DEGs
# Note: here we count DEG recovery rate: those found in full data, NOT all DEGs here
sample_raio = seq(0.1, 0.9, 0.1)

# init saturation results 
saturation_union_recover_rate = c()
saturation_intersect_recover_rate = c()

for (raio in sample_raio) {
  print(raio)
  # random sample 
  set.seed(123)
  select_1 = sample(sample_grp_1, size = round(raio * length(sample_grp_1)))
  select_2 = sample(sample_grp_2, size = round(raio * length(sample_grp_2)))
  keep_data = c(select_1, select_2)
  
  # filter df and metadata
  temp_df_count = input_data$counts[ , keep_data]
  temp_df_meta = input_data$group[keep_data, ]
  
  # call DEG, note we use "new_log2fc=1" here because we scale up the DEGs form full data

  temp_deg_deseq2 <- run_DESeq2_for_df_and_meta(df_count=temp_df_count, df_meta=temp_df_meta, single_condition=deg_condition)
  temp_deg_deseq2 <- filter_deg(temp_deg_deseq2, log2fc=cutoff_log2fc)
  
  temp_deg_edger <- run_edgeR_for_df_and_meta(df_count=temp_df_count, df_meta=temp_df_meta, single_condition=deg_condition)
  temp_deg_edger <- filter_deg(temp_deg_edger, log2fc=cutoff_log2fc)
  
  # create a new list to store results and nest it into out6
  temp_list = list()
  temp_list$deg_deseq2 = temp_deg_deseq2
  temp_list$deg_edger = temp_deg_edger
  temp_list$union_deg = union(rownames(temp_deg_deseq2), rownames(temp_deg_edger))
  temp_list$intersect_deg = intersect(rownames(temp_deg_deseq2), rownames(temp_deg_edger))
  temp_list$gt_recover_by_union = intersect(temp_list$union_deg, out6$overlap_deg)
  temp_list$gt_recover_by_union_ratio = length(temp_list$gt_recover_by_union) / length(out6$overlap_deg)
  temp_list$gt_recover_by_intersect = intersect(temp_list$intersect_deg, out6$overlap_deg)
  temp_list$gt_recover_by_intersect_ratio = length(temp_list$gt_recover_by_intersect) / length(out6$overlap_deg)
  
  # add plot stats
  saturation_union_recover_rate = c(saturation_union_recover_rate, temp_list$gt_recover_by_union_ratio)
  saturation_intersect_recover_rate = c(saturation_intersect_recover_rate, temp_list$gt_recover_by_intersect_ratio)
  
  # add this list to out6
  out6[[paste0("subsample_result_before_GEDE_ratio_", raio)]] <- temp_list
  rm(temp_list)
}


# save and plot results
saturation_results = list()
saturation_results$ratio = sample_raio
saturation_results$union_recover = saturation_union_recover_rate
saturation_results$intersect_recover = saturation_intersect_recover_rate

out6$saturation_before_GEDE = saturation_results
rm(saturation_results)

# save results
save(out6, file="output_real_data_analysis.Rdata")

# generate line plots
df <- data.frame(
  ratio = out6$saturation_before_GEDE$ratio,
  union_recover = out6$saturation_before_GEDE$union_recover,
  intersect_recover = out6$saturation_before_GEDE$intersect_recover
)

df_long <- melt(df, id.vars = "ratio", variable.name = "Type", value.name = "Recovery")

# Create the plot
p <- ggplot(df_long, aes(x = ratio, y = Recovery, color = Type, group = Type)) +
  geom_line(size = 1) +         # Draw lines
  geom_point(size = 2) +        # Add points
  geom_text(aes(label = round(Recovery, 3)), vjust = -0.5, size = 4) +  # Text labels
  labs(title = "Saturation Analysis", x = "Subsample ratio", y = "Recovery ratio") +
  scale_color_manual(values = c("red", "blue")) +  # Customize line colors
  theme_minimal()

# Save the plot to a PDF
ggsave("update_fig_7_saturation_analysis_before_GEDE.pdf", plot = p, width = 7, height = 5)



### 3.1, full data DEG after GEDE
############ Get DEGs for the full data
rm(deg_deseq2, deg_limma, deg_edger)
# DESeq2
deg_deseq2 <- run_DESeq2_for_df_and_meta(df_count=input_data$gede_counts, df_meta=input_data$group, single_condition=deg_condition, pvalue=cutoff_pvalue, log2fc=cutoff_log2fc)

# limma needs a voom transformed data
deg_limma <- run_limma_for_df_and_meta(df_count=input_data$gede_counts, df_meta=input_data$group, single_condition=deg_condition, pvalue=cutoff_pvalue, log2fc=cutoff_log2fc)

# edgeR
deg_edger <- run_edgeR_for_df_and_meta(df_count=input_data$gede_counts, df_meta=input_data$group, single_condition=deg_condition, pvalue=cutoff_pvalue, log2fc=cutoff_log2fc)


# while the number is too high, let's increase log2fc to 1~2 and then check the overlap and make a van di plot
filtered_deseq2 <- filter_deg(deg_deseq2, log2fc= cutoff_log2fc)
print(dim(filtered_deseq2))
filtered_limma <- filter_deg(deg_limma, log2fc = cutoff_log2fc)
print(dim(filtered_limma))
filtered_edger <- filter_deg(deg_edger, log2fc = cutoff_log2fc)
print(dim(filtered_edger))


# regular Venn
compare_deg_non_proportional_figure(rownames(filtered_deseq2), rownames(filtered_limma), rownames(filtered_edger), pdf_file = "update_venn_after_GEDE.pdf")

# add to results
out6$deg_gede_records_deseq2 <- deg_deseq2
out6$deg_gede_records_limma <- deg_limma
out6$deg_gede_records_edger <- deg_edger
out6$deg_sig_gede_deseq2 <- filtered_deseq2
out6$deg_sig_gede_limma <- filtered_limma
out6$deg_sig_gede_edger <- filtered_edger

out6$union_gede_deg <- union(rownames(filtered_deseq2), rownames(filtered_edger))   # no limma here for consistency
out6$overlap_gede_deg <- intersect(rownames(filtered_deseq2), rownames(filtered_edger)) 

# gede DEG recover and overlap 
out6$gede_union_recover_ratio = length(intersect(out6$union_gede_deg, out6$overlap_deg)) / length(out6$overlap_deg)
out6$gede_overlap_recover_ratio = length(intersect(out6$overlap_gede_deg, out6$overlap_deg)) / length(out6$overlap_deg)
out6$gede_extra_genes = setdiff(out6$overlap_gede_deg, out6$overlap_deg)

save(out6, file="output_real_data_analysis.Rdata")



### 3.2, saturation analysis after GEDE
group_values = levels(as.factor(df_meta[[deg_condition]]))
# split group and do random sample
sample_grp_1 = rownames(df_meta[df_meta[deg_condition] == group_values[1], ])
sample_grp_2 = rownames(df_meta[df_meta[deg_condition] == group_values[2], ])

# loop DEGs by different sample rate, note that for the 200 vs 300 RNA-seq case, 20% can recover 87.5% DEGs
# Note: here we count DEG recovery rate: those found in full data, NOT all DEGs here
sample_raio = seq(0.1, 0.9, 0.1)

# init saturation results 
saturation_union_recover_rate = c()
saturation_intersect_recover_rate = c()

for (ratio in sample_raio) {
  print(ratio)
  # random sample 
  set.seed(123)
  select_1 = sample(sample_grp_1, size = round(ratio * length(sample_grp_1)))
  select_2 = sample(sample_grp_2, size = round(ratio * length(sample_grp_2)))
  keep_data = c(select_1, select_2)
  
  # filter df and metadata (here is the data after GEDE)
  temp_df_count = input_data$gede_counts[ , keep_data]
  temp_df_meta = input_data$group[keep_data, ]
  
  # call DEG, note we use "new_log2fc=1" here because we scale up the DEGs form full data

  temp_deg_deseq2 <- run_DESeq2_for_df_and_meta(df_count=temp_df_count, df_meta=temp_df_meta, single_condition=deg_condition, pvalue=cutoff_pvalue, log2fc= new_log2fc)
  temp_deg_deseq2 <- filter_deg(temp_deg_deseq2, log2fc= cutoff_log2fc)
  
  temp_deg_edger <- run_edgeR_for_df_and_meta(df_count=temp_df_count, df_meta=temp_df_meta, single_condition=deg_condition, pvalue=cutoff_pvalue, log2fc= new_log2fc)
  temp_deg_edger <- filter_deg(temp_deg_edger, log2fc= cutoff_log2fc)
  
  # create a new list to store results and nest it into out6
  temp_list = list()
  temp_list$deg_deseq2 = temp_deg_deseq2
  temp_list$deg_edger = temp_deg_edger
  temp_list$union_deg = union(rownames(temp_deg_deseq2), rownames(temp_deg_edger))
  temp_list$intersect_deg = intersect(rownames(temp_deg_deseq2), rownames(temp_deg_edger))
  temp_list$gt_recover_by_union = intersect(temp_list$union_deg, out6$overlap_deg)
  temp_list$gt_recover_by_union_ratio = length(temp_list$gt_recover_by_union) / length(out6$overlap_deg)
  temp_list$gt_recover_by_intersect = intersect(temp_list$intersect_deg, out6$overlap_deg)
  temp_list$gt_recover_by_intersect_ratio = length(temp_list$gt_recover_by_intersect) / length(out6$overlap_deg)
  
  # add plot stats
  saturation_union_recover_rate = c(saturation_union_recover_rate, temp_list$gt_recover_by_union_ratio)
  saturation_intersect_recover_rate = c(saturation_intersect_recover_rate, temp_list$gt_recover_by_intersect_ratio)
  
  # add this list to out6
  out6[[paste0("subsample_result_after_GEDE_ratio_", ratio)]] <- temp_list
  rm(temp_list)
}



# save and plot results
saturation_results = list()
saturation_results$ratio = sample_raio
saturation_results$union_recover = saturation_union_recover_rate
saturation_results$intersect_recover = saturation_intersect_recover_rate

out6$saturation_after_GEDE = saturation_results
rm(saturation_results)



# save results
save(out6, file="output_real_data_analysis.Rdata")


# generate line plots
df <- data.frame(
  ratio = out6$saturation_after_GEDE$ratio,
  union_recover = out6$saturation_after_GEDE$union_recover,
  intersect_recover = out6$saturation_after_GEDE$intersect_recover
)

df_long <- melt(df, id.vars = "ratio", variable.name = "Type", value.name = "Recovery")

# Create line plot
p <- ggplot(df_long, aes(x = ratio, y = Recovery, color = Type, group = Type)) +
  geom_line(size = 1) +  # Draw lines
  geom_point(size = 2) + # Add points for better visibility
  geom_text(aes(label = round(Recovery, 3)), vjust = -0.5, size = 4) + # Add text labels
  labs(title = "Saturation Analysis", x = "Subsample ratio", y = "Recovery ratio") +
  scale_color_manual(values = c("red", "blue")) +  # Customize colors
  theme_minimal()

# Save the plot to a PDF
ggsave("update_fig_8_saturation_analysis_after_GEDE.pdf", plot = p, width = 7, height = 5)




# let's plot 2 saturations together
df_sat_overlap <- data.frame(
  ratio = out6$saturation_after_GEDE$ratio,
  raw_data = out6$saturation_before_GEDE$intersect_recover,
  gede_data = out6$saturation_after_GEDE$intersect_recover
)
df_long_overlap <- melt(df_sat_overlap, id.vars = "ratio", variable.name = "Type", value.name = "Recovery")


df_sat_union <- data.frame(
  ratio = out6$saturation_after_GEDE$ratio,
  raw_data = out6$saturation_before_GEDE$union_recover,
  gede_data = out6$saturation_after_GEDE$union_recover
)
df_long_union <- melt(df_sat_union, id.vars = "ratio", variable.name = "Type", value.name = "Recovery")

# plot them together
p1 <- ggplot(df_long_overlap, aes(x = ratio, y = Recovery, color = Type, group = Type)) +
  geom_line(size = 1) +  # Draw lines
  geom_point(size = 2) + # Add points for better visibility
  geom_text(aes(label = round(Recovery, 3)), vjust = -0.5, size = 4) + # Add text labels
  labs(title = "Saturation Analysis", x = "Subsample ratio", y = "Recovery ratio by OVERLAP DEGs") +
  ylim(0, 1) +
  scale_color_manual(values = c("red", "blue")) +  # Customize colors
  theme_minimal()

p2 <- ggplot(df_long_union, aes(x = ratio, y = Recovery, color = Type, group = Type)) +
  geom_line(size = 1) +  # Draw lines
  geom_point(size = 2) + # Add points for better visibility
  geom_text(aes(label = round(Recovery, 3)), vjust = -0.5, size = 4) + # Add text labels
  labs(title = "Saturation Analysis", x = "Subsample ratio", y = "Recovery ratio by UNION DEGs") +
  ylim(0, 1) +
  scale_color_manual(values = c("red", "blue")) +  # Customize colors
  theme_minimal()


# Save the figure as a file
combined_plot <- p1 + p2
ggsave(filename = "update_fig_9_compare_saturation_analysis.pdf", plot = combined_plot, width = 15, height = 6, dpi = 300)




# check differences
### Why genes are not DEG in subsample
extra_DEG_in_raw <- setdiff(out6$overlap_deg, out6$union_gede_deg)
# 3904 / 3303 for w.o. winsoratization
# 3732 / 3303 for w. winsoratization

# count overlaps
length(intersect(out6$overlap_deg, out6$union_gede_deg))  #  2677/2707, approx. 2/3

# make a overlap venn plot
compare_deg_non_proportional_figure(out6$union_deg, out6$overlap_deg, out6$union_gede_deg, names = c("Raw Union DEG", "Raw Overlap DEG", "GEDE Overlap DEG"), pdf_file = "update_venn_DEG_overlaps.pdf")

temp_record_deseq2 <- out6$deg_gede_records_deseq2[rownames(out6$deg_gede_records_deseq2) %in% extra_DEG_in_raw, ]
dim(temp_record_deseq2)  # 1061/915

# print results out
sum(temp_record_deseq2$padj > 0.05)   # 197
sum(abs(temp_record_deseq2$log2FoldChange) < 0.5) # 1060
temp_abs_log2fc <- abs(temp_record_deseq2$log2FoldChange)

pdf("update_fig_10_abs_log2fc_for_those_missing_in_gede.pdf", width = 7, height = 5)
hist(temp_abs_log2fc, 
     breaks = 10,  # Number of bins
     col = "blue",  # Color of bars
     border = "black",  # Edge color
     main = "Histogram of abs log2fc", 
     xlab = "Ratio", 
     xlim = c(0,1),
     ylab = "Frequency")
dev.off()

# we also plot those extra genes in GEDE
extra_gene_in_gede <- setdiff(out6$union_gede_deg, out6$overlap_deg)
temp_record_deseq2 <- out6$deg_all_records_deseq2[rownames(out6$deg_all_records_deseq2) %in% extra_gene_in_gede, ]
dim(temp_record_deseq2)  # 626
sum(temp_record_deseq2$padj > 0.05)   # 231
sum(abs(temp_record_deseq2$log2FoldChange) < 0.5) # 456

temp_abs_log2fc <- abs(temp_record_deseq2$log2FoldChange)

pdf("update_fig_11_abs_log2fc_for_those_extra_in_gede.pdf", width = 7, height = 5)
hist(temp_abs_log2fc, 
     breaks = 10,  # Number of bins
     col = "blue",  # Color of bars
     border = "black",  # Edge color
     main = "Histogram of abs log2fc", 
     xlab = "Ratio", 
     xlim = c(0,1),
     ylab = "Frequency")
dev.off()


# let's plot out LFC for those missing in GEDE
gede_non_gene_ids <- setdiff(out6$overlap_deg, out6$union_gede_deg)
length(gede_non_gene_ids)  # 1227 / 1025
gede_extra_gene_ids <- setdiff(out6$union_gede_deg, out6$overlap_deg)
length(gede_extra_gene_ids) # 626 / 596

log_df = log2(df_count+1)
dim(log_df)
dim(new_df)

# LFC (log-fold-changes) before/after GEDE, using log-data
n0 <- sum(numeric_group1==0) #191 black subjects
n1 <- sum(numeric_group1)    #298 white subjects
x <- ifelse(numeric_group1==0, -1/n0, 1/n1)
LFC.orig <- drop(log_df%*%x); LFC.gede <- drop(new_df%*%x)

pdf("update_fig_12_LFC_colored_special_points.pdf", width = 7, height = 7)
xy <- range(c(LFC.orig, LFC.gede))
# Plot all points in default color (black)
plot(LFC.gede ~ LFC.orig, xlim = xy, ylim = xy, 
     main = "LFC", xlab = "Original LFC", ylab = "GEDE-corrected LFC", pch = 20, col = "black")

# Add identity line
abline(0, 1, lty = 2)

# Now overplot specific points
# Assume names(LFC.orig) and names(LFC.gede) are gene IDs

# 1. Plot gede_non_gene_ids in red
common_red <- intersect(names(LFC.orig), gede_non_gene_ids)
points(LFC.orig[common_red], LFC.gede[common_red], pch = 20, col = "red")

# 2. Plot gede_extra_gene_ids in blue
common_blue <- intersect(names(LFC.orig), gede_extra_gene_ids)
points(LFC.orig[common_blue], LFC.gede[common_blue], pch = 20, col = "blue")

# Optionally add a legend
legend("topleft", legend = c("gede_non_gene_ids", "gede_extra_gene_ids"), 
       col = c("red", "blue"), pch = 20, bty = "n")

# Close PDF device
dev.off()





































































































