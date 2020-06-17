#!/usr/bin/env Rscript

## Compute TIE/power

##  0. Load libraries, parse arguments

library(optparse)
library(data.table)

option_list = list(
  make_option(c("-t", "--tests"), type="character", default=NULL,
              help="Testing results", metavar="character"),
  make_option(c("-l", "--level"), type="numeric", default=0.05,
              help="Significance level [default %default]", metavar="numeric"),
  make_option(c("-i", "--id_file"), type="character", default=NULL,
              help="IDs of variants generated under H1, if any", metavar="character"),
  make_option(c("-m", "--mtc"), type="character", default="none",
              help="Multiple testing correction: either 'none' or methods available in p.adjust [default %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output (MLM p-values) file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null (opt$tests) || is.null (opt$id_file) || is.null(opt$output) ){
  print_help(opt_parser)
  stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(123)

## 1. Load inputs
 tests <- fread(opt$tests, data.table = F)
 sel <- tryCatch({read.table(opt$id_file, as.is = T)[,1]}, 
                    error = function(e){NULL})

## 2. Filter and compute TIE/power. 
 if(!is.null(sel)){
   ids <- which(tests$rs %in% sel)
 } else {
   ids <- 1:nrow(tests) 
 }
 pv <- which(colnames(tests) %in% c("p_wald", "p_value"))

 if(opt$mtc != 'none'){ # Optional multiple testing correction
     tests[, pv] <- p.adjust(tests[, pv], method = opt$mtc)
 } 
 tie <- sum(tests[ids, pv] < opt$level) / length(ids) 

## 3. Store result
 write.table(tie, file = opt$output, col.names = F, row.names = F)
