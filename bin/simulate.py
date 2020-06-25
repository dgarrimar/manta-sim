#!/usr/bin/env python3

import gzip
import sys
import pickle
from random import randint, seed
from argparse import ArgumentParser

# Argument parsing
parser = ArgumentParser(description='Simulate genotypes with population stratification and relatedness.')
parser.add_argument("-g", "--genotypes", type=str, help="Genotypes in VCF from 1000g")
parser.add_argument("-A", "--ancestors", type=int, help="Number of ancestors")
parser.add_argument("-p", "--ancestors_pickle", type=str, help="Ancestor pickle file from 'generate_ancestors.py'")
parser.add_argument("-n", "--individuals", type=int, help="Number of individuals")
parser.add_argument("-b", "--blocksize", type=int, help="Number of variants per block")
parser.add_argument("-s", "--seed", type=int, help="Seed for random processes")
parser.add_argument("-o", "--out", type=str, help="Output VCF")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Rename parameters
A = args.ancestors
n = args.individuals
blocksize = args.blocksize

# Read pickled ancestors dict
pickle_in = open(args.ancestors_pickle,"rb")
ancestors = pickle.load(pickle_in)
pickle_in.close()

# Open input VCF
if args.genotypes.endswith('.gz'):
    fh = gzip.open(args.genotypes, 'rt')
else:
    fh = open(args.genotypes, 'r')

# Open output VCF
out_fh = open(args.out, 'w')

seed(args.seed)

# Read line by line
for line in fh:
    if line.startswith('##'):
        # Skip meta-information lines
        continue
    elif line.startswith('#CHROM'):
        # Read header
        header = line.strip().split("\t")
        count = 0 
    else:
        # For every block, select one ancestor per individual
        if (count % blocksize == 0):
            ancestors_block_idx = [header.index(ancestors[i][randint(0, A-1)]) for i in range(n)]
        # Obtain genotypes corresponding to the selected ancestor per individual
        line = line.strip().split('\t')
        gt = [line[idx] for idx in ancestors_block_idx]
        # Rename variants 
        line[2] = line[0] + "_" + line[1] + "_" + line[3] + "_" + line[4]
        # Remove INFO 
        line[7] = "."
        out_fh.write('\t'.join(line[0:9]) + '\t' + '\t'.join(gt) + '\n')
        count += 1

# Close files
fh.close()
out_fh.close()
