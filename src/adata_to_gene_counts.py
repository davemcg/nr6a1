import scanpy as sc
import pandas as pd
import argparse

adata = sc.read_h5ad(adata_path)
obs = adata.obs
obs['NR6A1_transform'] = adata[:,'ENSG00000148200'].to_df()

# move transformed counts into a layer
adata.layers['data'] = adata.X
# put the raw counts into X
adata.X = adata.raw.X
obs['NR6A1_raw'] = adata[:,'ENSG00000148200'].to_df()

obs.to_csv(out_file)
