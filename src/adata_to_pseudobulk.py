# on macbook use the scvi1.0.4 mamba env

import scanpy as sc
import pandas as pd
import argparse
from adpbulk import ADPBulk

parser = argparse.ArgumentParser(description = 'Use adp bulk to make pseudobulk matrix')
parser.add_argument("h5ad")
parser.add_argument("category", help = "comma separated (no whitespace) list of columns to aggregate against")
parser.add_argument("out_table")
args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad) 

the_cats = args.category.split(',')

adpb = ADPBulk(adata, the_cats, use_raw = True)
pseudobulk_matrix = adpb.fit_transform()
pseudobulk_matrix.to_csv(args.out_table)

