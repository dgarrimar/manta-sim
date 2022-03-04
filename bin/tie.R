#!/usr/bin/env Rscript

## Compute TIE/power

## 0. Load libraries, parse arguments

library(optparse)
library(data.table)

option_list = list(
    make_option(c("-t", "--tests"), type="character", default=NULL,
                help="Testing results", metavar="character"),
    make_option(c("-l", "--level"), type="numeric", default=0.05,
                help="Significance level [default %default]", metavar="numeric"),
    make_option(c("-m", "--mtc"), type="character", default="none",
                help="Multiple testing correction: either 'none' or methods available in p.adjust [default %default]", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output (adjusted p-values) file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null (opt$tests) || is.null(opt$output) ){
    print_help(opt_parser)
    stop("Required I/O files must be supplied\n", call.=FALSE)
}

set.seed(123)

## 1. Load inputs
tests <- fread(opt$tests, data.table = F)
tests <- tests[!is.na(tests[,2]), ] # rm NA's

## 2. Compute TIE/power. 

if(opt$mtc != 'none'){ # Optional multiple testing correction
    tests[, 2] <- p.adjust(tests[, 2], method = opt$mtc)
} 
tie <- mean(tests[, 2] < opt$level)  

## 3. Store result
write.table(tie, file = opt$output, col.names = F, row.names = F)
