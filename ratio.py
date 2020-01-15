#!/usr/bin/env python3

# Quick script to evaluate mat/pat ratio

import pandas as pd
import sys

MAXMUT = 5

mutfile = sys.argv[1]

data = pd.read_csv(mutfile, delim_whitespace=True, skipfooter=1) # skip last line bc my awk statement took too long to finish so file is truncated
#data = data[data["family"]==13976]
data = data[data["family"].apply(lambda x: x not in [11022, 12970, 13006, 13949, 13976, 14655, 11042, 11128, 11181, 13191, 13323, 14395])]

# Get rid of mutations occurring twice in same family
fam = data.groupby(["chrom","pos","family"], as_index=False).agg({"child": len})
fam = fam[fam["child"]==1].drop_duplicates()
data = pd.merge(data, fam[["chrom","pos","family"]], on=["chrom","pos","family"])

# Restrict to ctrls only
data = data[data["phenotype"]==1]

# Get rid of super mutable loci
loc = data.groupby(["chrom","pos"], as_index=False).agg({"child": len})
loc = loc[loc["child"]<MAXMUT]
data = pd.merge(loc[["chrom","pos"]], data, on=["chrom","pos"])

# Get rid of ones where we don't know POO
data = data[data["poocase"] != 4]

# Get rid of homozygous child genotypes
data["is.hom"] = data.apply(lambda x: len(set(x["child_gt"].split(",")))==1, 1)
data = data[~data["is.hom"]]

# Get ratio for each period
for period in range(1, 7):
    print("period=%s"%period)
    num_father = data[(data["period"]==period) & (data["poocase"]==2)].shape[0]
    num_mother = data[(data["period"]==period) & (data["poocase"]==3)].shape[0]
    print(num_father)
    print(num_mother)
    print(num_father/num_mother)

# Debug some
print(data[(data["poocase"]==2) & (data["period"]==2)][["chrom","pos","family","child","child_gt","mat_gt","pat_gt","poocase"]])
print(data[(data["poocase"]==3) & (data["period"]==2)][["chrom","pos","family","child","child_gt","mat_gt","pat_gt","poocase"]])
