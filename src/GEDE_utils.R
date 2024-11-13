# Stores all helper functions here for easier management
# Format: function + description + source from GEDE folder
print("GEDE_utils sourced!")
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
