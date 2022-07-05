```python
#展示celltype umap
#source activate python3
#python
import scdrs

import pandas as pd
import scanpy as sc

from anndata import AnnData

import matplotlib.pyplot as plt
import os

CORTEX_URL = "/share2/pub/zhenggw/zhenggw/scDRS/expression_mRNA_17-Aug-2014.txt"

# meta information
df_meta = pd.read_csv(CORTEX_URL, nrows=10, sep="\t", header=None)
df_meta = df_meta.iloc[:, 1:].T
columns = df_meta.iloc[0, :]
df_meta = df_meta.iloc[1:, :]
df_meta.columns = columns
df_meta = df_meta.set_index("cell_id")
df_meta.columns.name = None
df_meta["total mRNA mol"] = df_meta["total mRNA mol"].astype(float)

# expression information
df_expr = pd.read_csv(CORTEX_URL, skiprows=11, sep="\t", header=None).set_index(0)
df_expr.index.name = "gene"
# 1st column in backspin cluster, we don't need it here
df_expr = df_expr.iloc[:, 1:]
df_expr.columns = df_meta.index
df_expr = df_expr.T
raw_adata = AnnData(df_expr, obs=df_meta)
# raw_adata
# AnnData object with n_obs × n_vars = 3005 × 19972
# obs: 'tissue', 'group #', 'total mRNA mol', 'well', 'sex', 'age', 'diameter', 'level1class', 'level2class'

# assemble AnnData
sc.pp.filter_cells(raw_adata, min_genes=0)
sc.pp.filter_genes(raw_adata, min_cells=30)

# raw_adata
# AnnData object with n_obs × n_vars = 3005 × 14538
# obs: 'tissue', 'group #', 'total mRNA mol', 'well', 'sex', 'age', 'diameter', 'level1class', 'level2class', 'n_genes'
# var: 'n_cells'

raw_adata.raw = raw_adata

sc.pp.normalize_total(raw_adata, target_sum=1e4)
sc.pp.log1p(raw_adata)
sc.pp.highly_variable_genes(raw_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

raw_adata = raw_adata[:, raw_adata.var.highly_variable]

sc.pp.scale(raw_adata, max_value=10)
sc.tl.pca(raw_adata, svd_solver="arpack")
sc.pp.neighbors(raw_adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(raw_adata)
sc.tl.leiden(raw_adata)

#sc.pl.umap(raw_adata,color='tissue')

adata = raw_adata.raw.to_adata()
adata.obsp = raw_adata.obsp
adata.X = adata.X

# compile covariates
df_cov = pd.DataFrame(index=adata.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = adata.obs["n_genes"]

# load gene sets
df_gs = pd.read_excel(SUPP_TABLE_URL, sheet_name="ST6 MAGMA gene sets", index_col=0)
df_gs = df_gs.loc[["PASS_Schizophrenia_Pardinas2018", "UKB_460K.body_HEIGHTz"], :]
df_gs = scdrs.util.convert_gs_species(df_gs)
df_gs = df_gs.rename(
    index={"PASS_Schizophrenia_Pardinas2018": "SCZ", "UKB_460K.body_HEIGHTz": "Height"}
)

# save to files
if not os.path.exists("data/"):
    os.makedirs("data/")

adata.write_h5ad("data/expr.h5ad")
df_cov.to_csv("data/cov.tsv", sep="\t")
df_gs.reset_index().to_csv("data/geneset.gs", sep="\t", index=False)

```

```python
/share2/pub/zhenggw/zhenggw/anaconda3/envs/python3/bin/python /share2/pub/zhenggw/zhenggw/scDRS/compute_score.py \
    --h5ad_file /share2/pub/zhenggw/zhenggw/project1/test/data/expr.h5ad \
    --h5ad_species human \
    --gs_file /share2/pub/zhenggw/zhenggw/project1/test/data/geneset.gs \
    --gs_species human \
    --cov_file /share2/pub/zhenggw/zhenggw/project1/test/data/cov.tsv \
    --flag_filter True \
    --flag_raw_count True \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --out_folder /share2/pub/zhenggw/zhenggw/project1/test/data/
    
#####    
os.chdir('/share2/pub/zhenggw/zhenggw/project1/intest/')
adata = sc.read_h5ad("data/expr.h5ad")
df_gs = pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/umap/Bac_part.gs", sep="\t", index_col=0)

dict_score = {
trait: pd.read_csv(f"umap/{trait}.full_score.gz", sep="\t", index_col=0)
for trait in df_gs.index
}

for trait in dict_score:
adata.obs[trait] = dict_score[trait]["norm_score"]

sc.set_figure_params(figsize=[2.5, 2.5], dpi=150)
sc.pl.umap(
adata,
color="category",
ncols=1,
color_map="RdBu_r",
vmin=0,
vmax=4,
)

sc.pl.umap(
adata,
color=dict_score.keys(),
color_map="RdBu_r",
vmin=0,
vmax=4,
s=5,
save='.png',
)
```

