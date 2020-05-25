#!/usr/bin/env python

import pandas as pd
import gzip
import sys
from random import randint, seed
from argparse import ArgumentParser

# Argument parsing
parser = ArgumentParser(description='Simulate genotypes with population stratification and relatedness.')
parser.add_argument("-g", "--genotypes", type=str, help="Genotypes in VCF from 1000g")
parser.add_argument("-m", "--metadata", type=str, help="Population info from 1000g")
parser.add_argument("-A", "--ancestors", type=int, help="Number of ancestors")
parser.add_argument("-n", "--individuals", type=int, help="Number of individuals")
parser.add_argument("-b", "--blocksize", type=int, help="Number of variants per block")
parser.add_argument("-s", "--seed", type=int, help="Seed for random processes")
parser.add_argument("-o", "--outfile", type=str, help="Output also VCF format")
parser.add_argument("-v", "--outvcf", type=str, default=None, help="Output also VCF format")
parser.add_argument("-M", "--mode", type=str, help="Simulation mode. Only europeans: e, pop. stratification: p, both pe")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Read metadata
df = pd.read_csv(args.metadata, sep = '\t')
df = pd.DataFrame(df)

# Remove NaN columns
df = df.dropna(axis = 1)

# Define simulation parameters (simUnrelated,simRelated)
A = args.ancestors
n = args.individuals
blocksize = args.blocksize

ancestors = dict()

# Define simulation mode
if 'e' in args.mode:
    europeans = True
else:
    europeans = False
if 'p' in args.mode:
    popstr = True
else:
    popstr = False

seed(0)

if europeans:
    if popstr:
        pops = list(df.loc[df['super_pop']=='EUR', 'pop'].unique())
        for i in range(n):
            pop = pops[randint(0, len(pops)-1)]
            inds = list(df.loc[df['pop'] == pop, 'sample'])
            ancestors[i] = [inds[randint(0, len(inds)-1)] for i in range(A)]
          
    else:
        inds = list(df.loc[df['super_pop'] == 'EUR', 'sample'])
        for i in range(n):
            ancestors[i] = [inds[randint(0, len(inds)-1)] for i in range(A)]
else:
    if popstr:
        pops = list(df['pop'].unique())
        for i in range(n):
            pop = pops[randint(0, len(pops)-1)]
            inds = list(df.loc[df['pop'] == pop, 'sample'])
            ancestors[i] = [inds[randint(0, len(inds)-1)] for i in range(A)]
    else:
        inds = list(df['sample'])
        for i in range(n):
            ancestors[i] = [inds[randint(0, len(inds)-1)] for i in range(A)]

# Format output
def format_out(line, genotypes, out_fh, outvcf_fh):
    for i, gt in enumerate(genotypes):
        if (gt == '0|0'):
            genotypes[i] = "0"
        elif(gt == '1|0' or genotypes[i] == '0|1'):
            genotypes[i] = "1"
        elif(gt == '1|1'):
            genotypes[i] = "2"
    if (line[2] == "."):
        snp = line[0] + "_" + line[1]
    else:
        snp = line[2]
    out_fh.write(snp + ", " + ", ".join(line[3:5]) + ", " + ", ".join(genotypes) + '\n')
    if outvcf_fh is not None:
        outvcf_fh.write('\t'.join(line[0:9]) + '\t' + '\t'.join(genotypes) + '\n')
        
# Open VCF
if args.genotypes.endswith('.gz'):
    fh = gzip.open(args.genotypes, 'rt')
else:
    fh = open(args.genotypes, 'r')

# Open output file(s)
out_fh = open(args.outfile, 'w')
if args.outvcf is not None:
    outvcf_fh = open(args.outvcf, 'w')
else:
    outvcf_fh = None

seed(args.seed)

# Read line by line
for line in fh:
    if line.startswith('##'):
        # skip comment lines
        continue
    elif line.startswith('#CHROM'):
        # read header
        header = line.strip().split("\t")
        count = 0 
    else:
        # For every block, select one ancestor per individual
        if (count % blocksize == 0):
            ancestors_block_idx = [header.index(ancestors[i][randint(0, A-1)]) for i in range(n)]
        # Obtain genotypes corresponding to the selected ancestor per individual
        line = line.strip().split('\t')
        gt = [line[idx] for idx in ancestors_block_idx]
        format_out(line, gt, out_fh, outvcf_fh)
        count += 1

# Close files
fh.close()
out_fh.close()
if(args.outvcf is not None):
    outvcf_fh.close()
