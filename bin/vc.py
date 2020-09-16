#!/usr/bin/env python

## Import modules
import sys
from argparse import ArgumentParser
from pylmm import input
from pylmm.lmm import LMM
from scipy import linalg
import numpy as np
import time

## Define options
def define_options():
    # Argument parsing
    parser = ArgumentParser(description = 'Calculate variance components')
    parser.add_argument("-b", "--bfile", type = str,
        help = "PLINK binary bed file prefix")
    parser.add_argument("-p", "--pheno", type = str,
        help = "Phenotype file in PLINK format: FID, IID, [PHENOTYPES]")
    parser.add_argument("-k", "--kinship", type = str,
        help = "Kinship matrix (n x n plain text file)")
    parser.add_argument("-c", "--covfile", type = str,
        help = "Covariate file in PLINK format: FID, IID, [COVARIATES]")
    parser.add_argument("-o", "--outfile", type = str,
        help = "Output file")
    parser.add_argument("-v", "--verbose", action = "store_true",
        help = "Print extra info [default=%(default)s]", default = False)
        ## Other options
        # remove missing genotypes
        # no mean
    return parser

parser = define_options()
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

## Read PLINK input (and covariates)

IN = input.plink(args.bfile, type = 'b', phenoFile = args.pheno)

if args.covfile:
    C = IN.getCovariates(args.covfile) # Read covariate file, write into input.plink
    X0 = np.hstack([np.ones((IN.phenos.shape[0], 1)), C]) # Added global mean
    if np.isnan(X0).sum():
        parser.error("The covariate file should not contain missing values.")
    if args.verbose: sys.stderr.write("Read %d covariate(s) from the covariate file\n" % C.shape[1])
else:
    X0 = np.ones((IN.phenos.shape[0], 1))

if np.isnan(IN.phenos).sum():
    parser.error("The phenotype file should not contain missing values.")

## Read Kinship: this method seems to be the fastest and works if you already know the size of the matrix

begin = time.time()
if args.kinship[-3:] == '.gz':
    import gzip
    f = gzip.open(args.kinship,'r')
    F = f.read() # Might exhaust mem if the file is very large
    K = np.fromstring(F, sep = ' ') # Assume space separated
    f.close()
else:
    K = np.fromfile(open(args.kinship, 'r'), sep = ' ')
    K.resize((len(IN.indivs), len(IN.indivs)))
end = time.time()
if args.verbose: sys.stderr.write("Read the %d x %d kinship matrix in %0.3fs\n" % (K.shape[0], K.shape[1], end-begin))

## Compute variance components
phenoNum = IN.phenos.shape[1]
sys.stderr.write("Read %d phenotype(s)\n" % phenoNum)
out = open(args.outfile, 'w')
begin = time.time()
Kva, Kve = linalg.eigh(K)
end = time.time()
if args.verbose: sys.stderr.write("Kinship eigendecomposition time: %0.3fs\n" % (end - begin))
for i in range(0, phenoNum):
    if args.verbose: sys.stderr.write("Phenotype number: %d \n" % (i+1))
    Y = IN.phenos[:,i]
    L = LMM(Y, K, Kva, Kve, X0, verbose = args.verbose)
    if args.verbose: sys.stderr.write("Computing fit for null model\n")
    L.fit()
    if args.verbose:
        sys.stderr.write("\t heritability=%0.5f, varG=%0.5f, varE=%0.5f\n" % (L.optH,L.optH*L.optSigma, L.optSigma*(1-L.optH)))
        out.write("%0.5f\t%0.5f\n" % (L.optH*L.optSigma, L.optSigma*(1-L.optH)))
out.close()
