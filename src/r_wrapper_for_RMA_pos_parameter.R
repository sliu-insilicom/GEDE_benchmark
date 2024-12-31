date()
library(rmarkdown)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No Rmd file provided. Usage: Rscript render_script.R your_file.Rmd")
}
rmd_file <- args[1]

rmarkdown::render(rmd_file, output_format = "html_document")

date()
