#!/usr/bin/env python

import pandas as pd
import sys
import pickle
from random import randint, seed
from argparse import ArgumentParser

# Argument parsing
parser = ArgumentParser(description='Generate ancestors for subsequent simulation with simulate.py.')
parser.add_argument("-m", "--metadata", type=str, help="Population info from 1000g")
parser.add_argument("-A", "--ancestors", type=int, help="Number of ancestors")
parser.add_argument("-n", "--individuals", type=int, help="Number of individuals")
parser.add_argument("-s", "--seed", type=int, help="Seed for random processes")
parser.add_argument("-p", "--pickle", type=str, help="Ancestors pickle file for 'simulate.py'")
parser.add_argument("-M", "--mode", type=str, help="Simulation mode. Only europeans: e, pop. stratification: p, both: pe")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Read metadata
df = pd.read_csv(args.metadata, sep = '\t')
df = pd.DataFrame(df).dropna(axis = 1)

# Rename parameters
A = args.ancestors
n = args.individuals

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

seed(args.seed)

if europeans:
    if popstr:
        pops = list(df.loc[df['super_pop']=='EUR', 'pop'].unique())
        for i in range(n):
            pop = pops[randint(0, len(pops)-1)]
            inds = list(df.loc[df['pop'] == pop, 'sample'])
            ancestors[i] = [inds[randint(0, len(inds)-1)] for j in range(A)]
          
    else:
        inds = list(df.loc[df['super_pop'] == 'EUR', 'sample'])
        for i in range(n):
            ancestors[i] = [inds[randint(0, len(inds)-1)] for j in range(A)]
else:
    if popstr:
        pops = list(df['pop'].unique())
        for i in range(n):
            pop = pops[randint(0, len(pops)-1)]
            inds = list(df.loc[df['pop'] == pop, 'sample'])
            ancestors[i] = [inds[randint(0, len(inds)-1)] for j in range(A)]
    else:
        inds = list(df['sample'])
        for i in range(n):
            ancestors[i] = [inds[randint(0, len(inds)-1)] for j in range(A)]

# Pickled output file
pickle_out = open(args.pickle,"wb")
pickle.dump(ancestors, pickle_out)
pickle_out.close()
