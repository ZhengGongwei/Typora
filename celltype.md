```python
#Intest

#task：
	#保存稀疏矩阵 h5ad文件 后续scDRS输入文件
	#cell tpye umap图
    
import scanpy as sc
import pandas as pd
import os
import loompy
import matplotlib.pyplot as plt

raw_adata = sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
# assemble AnnData
sc.pp.filter_cells(raw_adata, min_genes=0)
sc.pp.filter_genes(raw_adata, min_cells=30)
raw_adata.raw = raw_adata

#sc.pl.umap(raw_adata,color='tissue')

sc.pp.normalize_total(raw_adata, target_sum=1e4)
sc.pp.log1p(raw_adata)
sc.pp.highly_variable_genes(raw_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

raw_adata = raw_adata[:, raw_adata.var.highly_variable]

sc.pp.scale(raw_adata, max_value=10)
sc.tl.pca(raw_adata, svd_solver="arpack")
sc.pp.neighbors(raw_adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(raw_adata)
sc.tl.leiden(raw_adata)

adata = raw_adata.raw.to_adata()
adata.obsp = raw_adata.obsp
adata.X = adata.X

raw_adata = sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
sc.pp.filter_cells(raw_adata, min_genes=0)
sc.pp.filter_genes(raw_adata, min_cells=30)
raw_adata.raw = raw_adata


/share2/pub/zhenggw/zhenggw/anaconda3/envs/python3/bin/python /share2/pub/zhenggw/zhenggw/scDRS/compute_score.py \
    --h5ad_file /share2/pub/zhenggw/zhenggw/project1/test/data/expr.h5ad \
    --h5ad_species mouse \
    --gs_file /share2/pub/zhenggw/zhenggw/project1/test/data/geneset.gs \
    --gs_species mouse \
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

