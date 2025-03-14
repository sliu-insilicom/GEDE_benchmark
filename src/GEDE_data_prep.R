# Stores all helper functions here for easier management
# Format: function + description + source from GEDE folder
print("GEDE_data_prep sourced!")
library(DESeq2)
library(EnhancedVolcano)
# library(SummarizedExperiment)

###### from file: Original_data/data_prep.R
# link: https://www.dropbox.com/home/new_GEDE/Core_functions/Original_Data?preview=data_prep.R


# This function reads in a Rdata file storing a RSE object 
# @parameter: 
#     path to Rdata file containing a RSE object, which MUST have name "rse_gene"
#     metaVar: a single metadata item used to assign group variables
#     to_remove: a single group label used to filter input data
# @output:
#     a list containing: 1) group variable, 2) scaled counts table

read_data <- function(path, metaVar = NULL, to_remove = NULL) {
  load(path)
  # old code:
  # rse_gene <- scale_counts(rse_gene)  # this seems added a new array into object, but in re3 directly return the array
  # counts <- assays(rse_gene)$counts
  # new:
  counts <- transform_counts(rse_gene)
  if (!is.null(metaVar)) {
    # suggest for multiple groups if using metaVar as a vector:
    # groups <- lapply(metaVars, function(var) eval(parse(text = paste0("rse_gene@colData@listData$", var))))
    group <- eval(parse(text = paste0("rse_gene@colData@listData$", metaVar)))
    # if want to use as a vector
    # to_remove <- c(which(is.na(group)), which(group %in% to_remove))
    to_remove <- c(which(is.na(group)), which(group == to_remove))
    if (length(to_remove) > 0) {
      group <- group[-to_remove]
      counts <- counts[,-to_remove]
    }
  } else {
    group <- rep("NA", ncol(counts))
  }
  rm(rse_gene)
  return(list(group = group, counts = counts))
}



# This fuction makes the metadata variable of interest binary for DEA
# @parameter:
#     Data: The output list object from function "read_data"
#     positive_name: one value in the group variable 
#     negative_name: given value to set all non-positive values
#     to_remove: used to subset data
# @output:
#     still the list object but binarize the group variable

metadata_binarize <- function(Data, positive_name, negative_name, to_remove = NULL) {
  # use to_remove indices to remove some samples
  if (!is.null(to_remove)) {
    Data$counts <- Data$counts[,-to_remove]
    Data$group <- Data$group[-to_remove]
  }
  Data$group[Data$group != positive_name] <- negative_name
  return(Data)
}



# This function preps the list object ready for DEA, group must be binary
# @parameter
#     Data: The list object from metadata_binarize output (binary group)
# @output
#     Same object but added a design component

dea_prep <- function(Data) {
  group <- Data$group
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c(gsub(" ", "", colnames(design)[1], fixed = TRUE), 
                        gsub(" ", "", colnames(design)[2], fixed = TRUE))
  Data$design <- design
  ## define contrast (binary)
  contrast <- makeContrasts(contrasts = paste(gsub(" ", "", colnames(design)[1], fixed = TRUE), 
                                              "-", 
                                              gsub(" ", "", colnames(design)[2], fixed = TRUE)), 
                            levels = design)
  Data$contrast <- contrast
  Data$group <- paste0("group", gsub(" ", "", Data$group))
  return(Data)
}



# This function subset gene records by a vector of gene ID 
# (SL: I simplified the input by removing filepath but just using a vector for subseting)
# @parameter
#     Data: The list object previous functions (not necessarily DEA ready)
#     zero_threshold: cutoffs to remove gene records if missing portion is higher than this
#     gene_vector: a vector of gene IDs used to subset input data
#     include_voom: (limma) Transform count data to log2 counts-per-million (logCPM)
# @output:
#     Same list object with filtering and add a new component of voom transformed matrix "voom_counts"

gene_selection <- function(Data, zero_threshold = 0.9, # 1 means do not filter
                           gene_vector = NULL, include_voom = TRUE) {
  if (!is.null(gene_vector)) {
    Data$counts <- Data$counts[gene_vector,]
  }
  if (zero_threshold < 1) {
    ## get the ratio of zeros of each gene
    ratio_zeros <- rowSums(Data$counts == 0)/ncol(Data$counts)
    ## remove genes with ratio_zeros > zero_threshold
    Data$counts <- Data$counts[ratio_zeros <= zero_threshold, ]
  }
  if (include_voom) {
    y <- voom(Data$counts, plot = FALSE)
    Data$voom_counts <- y$E
  }
  return(Data)
}




# Check normality of rows (gene expression arcoss samples) of an input df
check_normality <- function(Data, pvalue=0.05) {
  # remove rows with no variation
  filtered_data <- Data[apply(Data, 1, function(row) var(row) > 0), ]
  # test for normality
  shapiro_results <- apply(filtered_data, 1, function(row) shapiro.test(row)$p.value)
  # adjust for multiple test
  adjusted_pvalues <- p.adjust(shapiro_results, method = "fdr")
  # extract sig genes
  unnormal_genes <- adjusted_pvalues[which(adjusted_pvalues < pvalue)]
  return(unnormal_genes)
}




# pull 1 dataset from recount3 by its keyword, e.g. "LUAD", "BRCA", and extract group variables
# @parameter
# available_project: project list obtained by recount3 available_project, has a fixed format
# keyword: unique keyword to extract 1 database, e.g. "LUAD"
# group_vector: a vector of clinical variables to keep as group variables

pull_recount3_data <- function(available_project, keyword, group_vector) {
  print(keyword)
  print(group_vector)
  proj_rse = subset(selected_proj, project == keyword & project_home == "data_sources/tcga")
  rse_gene <- create_rse(proj_rse)
  # filter non-coding gene
  if ("gene_type" %in% names(mcols(rowRanges(rse_gene)))) {
    coding_genes <- rowRanges(rse_gene)$gene_type == "protein_coding"
    rse_gene <- rse_gene[coding_genes, ]
  }
  # dim
  print(dim(rse_gene))
  # extract df and group information
  counts = compute_read_counts(rse_gene)
  groups <- as.data.frame(colData(rse_gene)[, group_vector])
  # return a list object
  return(list(group = groups, counts = counts))
}




# run DESeq2 for a raw count matrix, metadata df, and metadata colname
# need a binary "single_condition"
run_DESeq2 <- function(df_count, df_meta, single_condition, pvalue=0.05, log2fc=0.5) {
  # build object
  dds <- DESeqDataSetFromMatrix(
    countData = df_count,
    colData = df_meta,
    design = as.formula(paste0("~ ", single_condition))  # Specify the group variable
  )
  # filter low counts (don't use for now because we didn't include it in edgeR or limma)
  # dds <- dds[rowSums(counts(dds)) > 1000, ]  # Keep genes with a total count >1000, as we have hundreds of samples
  # run DEG
  dds <- DESeq(dds)
  compare_conditions <- levels(as.factor(df_meta[[single_condition]]))
  res <- results(dds, contrast = c(single_condition, compare_conditions[1], compare_conditions[2]))
  print(summary(res))
  # Filter significant genes
  significant_genes  <- res[!is.na(res$padj) & res$padj < pvalue & res$log2FoldChange>log2fc, ]
  print(dim(significant_genes))
  # make some plots
  EnhancedVolcano(res,
                  lab = rownames(res),                   # Gene labels
                  x = 'log2FoldChange',                  # Fold change
                  y = 'pvalue',                          # p-value
                  title = 'Volcano Plot',
                  pCutoff = pvalue,                        # P-value cutoff
                  FCcutoff = log2fc                           # Fold change cutoff
  )
  plotMA(res, main = "MA Plot", ylim = c(-5, 5))
  # return results
  return(significant_genes)
  
}


# run limma for a voom transformed matrix, metadata vector (binary)
run_limma <- function(df_transformed, meta_vector, pvalue=0.05, log2fc=0.5, betahat = 0.5) {
  meta_vector <- factor(meta_vector)
  design <- model.matrix(~ meta_vector)
  fit <- lmFit(df_transformed, design)
  fit2 <- eBayes(fit)
  results <- topTable(fit2, number=Inf, sort.by="none", adjust.method="BH")
  # add betahat for simulation
  coefficients <- fit2$coefficients
  combined_results <- cbind(results, Coefficients=coefficients)
  # may adjust criteria later
  DEGs <- which(combined_results$adj.P.Val < pvalue & abs(combined_results$Coefficients.meta_vector1) > betahat)
  #selected_genes <- rownames(combined_results)[DEGs]
  return(DEGs)
}











###################################################################################################
# train/test/valid split
train_valid_test_split <- function(data_group, train_ratio = 0.4, include_valid = TRUE, valid_ratio = 2/3) {
  all_index <- seq(1,length(data_group))
  train_index <- createDataPartition(data_group, p=train_ratio, list=FALSE, times=1)
  if (include_valid) {
    valid_index <- createDataPartition(data_group[-train_index], p=valid_ratio, list=FALSE, times=1)
    valid_index <- all_index[-train_index][valid_index]
    test_index <- all_index[-c(train_index, valid_index)]
    return(list(train_index = train_index, valid_index = valid_index, test_index = test_index))
  } else {
    test_index <- all_index[-train_index]
    return(list(train_index = train_index, test_index = test_index))
  }
}


###################################################################################################
# small test sets
sample_small_test <- function(group_test, nEach = 20, nIter = 50) {
  stest_indices <- list()
  for (i in 1:nIter) {
    set.seed(i)
    stest_index <- createDataPartition(group_test, p=nEach/length(group_test), 
                                       list=FALSE, times=1)
    stest_indices[[i]] <- stest_index
  }
  return(stest_indices)
}


#stest_indices <- sample_small_test(Data$group[test_index])
#save(file = 'train_valid_test_index.RData', train_index, valid_index, test_index, stest_indices)



