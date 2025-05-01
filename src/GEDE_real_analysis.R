# store all helper functions used for the real data analysis
print("GEDE_real_analysis sourced!")

# consider merge all 3 Rscrip together later


#### DEG analysis: assume df_count row for gene + col for sample; df_meta may contain multiple columns but only use 1
# run DESeq2 (similarly defined in "GEDE_data_prep.R", here we use same input for three so reset)
run_DESeq2_for_df_and_meta <- function(df_count, df_meta, single_condition, pvalue=0.05, log2fc=0.5) {
  # build object
  dds <- DESeqDataSetFromMatrix(
    countData = df_count,
    colData = df_meta,
    design = as.formula(paste0("~ ", single_condition))  # Specify the group variable
  )
  # filter low counts, this is a subjective cutoff can determined arbitrarily
  drop <- which(apply(cpm(dds), 1, max) < 1)
  dds <- dds[-drop, ]
  # run DEG
  dds <- DESeq(dds)
  compare_conditions <- levels(as.factor(df_meta[[single_condition]]))
  res <- results(dds, contrast = c(single_condition, compare_conditions[1], compare_conditions[2]))
  # return the full table for DEG checking
  return(res)
  
  # print(summary(res))
  # Filter significant genes
  # significant_genes  <- res[!is.na(res$padj) & res$padj < pvalue & abs(res$log2FoldChange)>log2fc, ]
  # print(dim(significant_genes))
  # 
  # return(significant_genes)
}


# run_limma (similarly defined in "GEDE_data_prep.R", here we use same input for three so reset)
run_limma_for_df_and_meta <- function(df_count, df_meta, single_condition, pvalue=0.05, log2fc=0.5) {
  # setup design
  temp_metadata <- df_meta[colnames(df_count), , drop = TRUE]  # Reorder to match count data
  temp_group <- factor(temp_metadata[[single_condition]])
  levels(temp_group) <- gsub(" ", "_", levels(temp_group))  # replace space by _ for the limma contrast
  design <- model.matrix(~temp_group)
  # normalize data
  # Convert raw counts to DGEList object (for voom normalization)
  dge <- DGEList(counts = df_count)
  dge <- calcNormFactors(dge)  # Normalize for library size
  # filter low expressed genes (CPM <1 )
  drop <- which(apply(cpm(dge), 1, max) < 1)
  dge <- dge[-drop,] 
  voom_count <- voom(dge, design, plot = FALSE)
  # Fit a linear model
  fit <- lmFit(voom_count)
  fit2 <- eBayes(fit)

  # Extract differentially expressed genes (DEGs)
  res <- topTable(fit2, number=Inf, sort.by = "P")
  # return the full table for DEG checking
  return(res)
  # significant_genes  <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < pvalue & abs(res$logFC)>log2fc, ]
  
  # return(significant_genes)
}


# run edgeR
run_edgeR_for_df_and_meta <- function(df_count, df_meta, single_condition, pvalue=0.05, log2fc=0.5) {
  # setup design
  temp_metadata <- df_meta[colnames(df_count), , drop = TRUE]  # Reorder to match count data
  temp_group <- factor(temp_metadata[[single_condition]])
  levels(temp_group) <- gsub(" ", "_", levels(temp_group))  # replace space by _ for the limma contrast
  design <- model.matrix(~temp_group)
  # setup DGE object
  d <- DGEList(counts=df_count,group=temp_group)
  # filter low expressed data
  drop <- which(apply(cpm(d), 1, max) < 1)
  d <- d[-drop,] 
  # normalization
  d <- calcNormFactors(d)
  # est dispersion
  d1 <- estimateDisp(d, design)
  # Fit GLM (generalized linear model)
  fit <- glmFit(d1, design)
  # Define contrast (extracting second coefficient)
  lrt <- glmLRT(fit, coef = 2)
  
  # Extract results table
  res <- topTags(lrt, n = Inf)$table  # Get all genes
  # return the full table for DEG checking
  return(res)
  
  # Filter significant DEGs (adjusted p-value < pvalue & absolute log2FC > log2fc)
  # significant_genes  <- res[!is.na(res$FDR) & res$FDR < pvalue & abs(res$logFC)>log2fc, ]
  # 
  # return(significant_genes)
}


# filter deg list by differnt padj and log2fc cutoff
filter_deg <- function(input_deg, pvalue=0.05, log2fc=1) {
  # determine input data from DESeq2/limma/edgeR
  res = input_deg # just for convenience
  if ("padj" %in% colnames(input_deg)) {   # DESeq2
    out <- res[!is.na(res$padj) & res$padj < pvalue & abs(res$log2FoldChange)>log2fc, ]
  }
  else if ("adj.P.Val" %in% colnames(input_deg)) {   # limma
    out <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < pvalue & abs(res$logFC)>log2fc, ]
  }
  else if ("FDR" %in% colnames(input_deg)) {  # edgeR
    out <- res[!is.na(res$FDR) & res$FDR < pvalue & abs(res$logFC)>log2fc, ]
  }
  else {
    print("Unexpected data type!!!")
    return()
  }
  
  return(out)
}



# check overlap and make venn diagram
compare_deg_non_proportional_figure <- function(v1, v2, v3, names=c("DESeq2", "limma", "EdgeR")) {
  # proportion-irrelavant figure
  # Compute overlaps
  n12 <- length(intersect(v1, v2))  # v1 ∩ v2
  n13 <- length(intersect(v1, v3))  # v1 ∩ v3
  n23 <- length(intersect(v2, v3))  # v2 ∩ v3
  n123 <- length(Reduce(intersect, list(v1, v2, v3)))  # v1 ∩ v2 ∩ v3
  
  # Generate the Venn diagram
  venn.plot <- draw.triple.venn(
    area1 = length(v1),
    area2 = length(v2),
    area3 = length(v3),
    n12 = n12, 
    n13 = n13, 
    n23 = n23, 
    n123 = n123, 
    category = names,  # Assign method names
    fill = c("red", "blue", "green"),  # Custom colors
    alpha = 0.5, 
    lwd = 2, 
    cex = 1.5, 
    cat.cex = 1.5, 
    cat.col = c("red", "blue", "green")
  )
  
  grid.newpage()  # Create a new plotting page
  grid.draw(venn.plot)  # Draw the Venn diagram
}


compare_deg_proportional_figure <- function(v1, v2, v3, names=c("DESeq2", "limma", "EdgeR")) {
  # venn based on proportion
  venn_data <- list(
    "DESeq2" = v1,
    "limma" = v2,
    "EdgeR" = v3
  )
  
  fit <- euler(venn_data)
  
  p <- plot(fit, 
         fills = c("red", "blue", "green"), 
         edges = TRUE, 
         labels = TRUE,
         quantities = TRUE,  # Show values
         main = "Proportional Venn Diagram")
  
  print(p)
  return(fit)
}


