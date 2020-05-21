#!/usr/bin/env python

import pandas as pd
import gzip
import sys
from random import randint, seed
from argparse import ArgumentParser

# Argument parsing
parser = ArgumentParser(description='Simulate genotypes with population stratification and relatedness.')
parser.add_argument("-g", "--genotypes", type=str, help="Genotypes in VCF from 1000g")
parser.add_argument("-p", "--popdata", type=str, help="Population info from 1000g")
parser.add_argument("-A", "--ancestors", type=int, help="Number of ancestors")
parser.add_argument("-n", "--individuals", type=int, help="Number of individuals")
parser.add_argument("-b", "--blocksize", type=int, help="Number of variants per block")
parser.add_argument("-s", "--seed", type=int, help="Seed for random processes")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Read metadata
df = pd.read_csv(args.popdata, sep = '\t')
df = pd.DataFrame(df)

# Remove NaN columns
df = df.dropna(axis = 1)

# Define simulation parameters (simUnrelated,simRelated)
A = args.ancestors
n = args.individuals
blocksize = args.blocksize
EUR = list(df.loc[df['super_pop']=='EUR', 'sample'])
POP = EUR

if args.genotypes.endswith('.gz'):
    fh = gzip.open(args.genotypes, 'rt')
else:
    fh = open(args.genotypes, 'r')

seed(args.seed)
for line in fh:
    if line.startswith('##'):
        # skip comment lines
        continue
    elif line.startswith('#CHROM'):
        # read header, generate ancestor list per individual
        header = line.strip().split("\t")
        ancestors = dict()
        for i in range(n):
            ancestors[i] = [POP[randint(0, len(POP)-1)] for i in range(A)]
        count = 0
        # print ('\t'.join(header[0:9]) + '\t' + '\t'.join(['sim' + str(i+1) for i in range(n)]))
    else:
        # For every block variants, change ancestors
        if (count % blocksize == 0):
            ancestors_block_idx = [header.index(ancestors[i][randint(0, A-1)]) for i in range(n)]
        line = line.strip().split('\t')
        gt = [line[idx] for idx in ancestors_block_idx]
        # Obtain genotypes corresponding to the selected ancestors
        print('\t'.join(line[0:9]) + '\t' + '\t'.join(gt))
        # print('\t'.join(gt))
        count += 1

fh.close()
