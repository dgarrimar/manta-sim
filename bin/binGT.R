#!/usr/bin/env Rscript

## Simulate binomial SNPs with a given MAF

## 0. Load libraries, parse arguments, define functions

library(optparse)

option_list = list(
    make_option(c("-s", "--seed"), type="numeric", default=0,
                help="Set seed for random processes [default %default]", metavar="numeric"),
    make_option(c("-n", "--nb_samples"), type="numeric",
                help="Total number of samples", metavar="numeric"),
    make_option(c("-m", "--maf"), type="numeric",
                help="MAF", metavar="numeric"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output (simulated genotypes, VCF format) file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$n) || is.null(opt$maf) || is.null(opt$output)){
    print_help(opt_parser)
    stop("Required I/O files must be supplied\n", call.=FALSE)
}

## 1. Simulate

set.seed(opt$seed)

X <- rbinom(opt$n, 2, opt$maf)
X[X==0] <- "0/0"
X[X==1] <- "0/1"
X[X==2] <- "1/1"

X <- c(1, 1, "dummy", "A", "T", ".", ".", "PR", "GT", X)
names(X) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
              "FILTER", "INFO", "FORMAT", paste0("S", 1:opt$n))

## 2. Save

write.table(t(X), file = opt$output, col.names = T, row.names = F, sep = "\t", quote = F)

