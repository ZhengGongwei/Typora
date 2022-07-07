import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=50, dpi_save=300, figsize=(5, 5))

# load filtered data 
# sc.pp.filter_cells(adata, min_genes=250) sc.pp.filter_genes(adata, min_cells=50)
adata=sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
adata=sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest.h5ad")
# AnnData object with n_obs × n_vars = 428469 × 24887

#add bac score info
#Epithelial
bac11=pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/data/Bac_11.genes.score.gz", sep="\t", index_col=0)
adata.obs["bac11"] = bac11["norm_score"]

#Epithelial (?Neuronal)
bac12=pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/data/Bac_12.genes.score.gz", sep="\t", index_col=0)
adata.obs["bac12"] = bac12["norm_score"]

#Epithelial Endothelial Mesenchymal
bac25=pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/data/Bac_25.genes.score.gz", sep="\t", index_col=0)
adata.obs["bac25"] = bac25["norm_score"]

#Epithelial
bac27=pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/data/Bac_27.genes.score.gz", sep="\t", index_col=0)
adata.obs["bac27"] = bac27["norm_score"]

#Red blood cells
bac36=pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/data/Bac_36.genes.score.gz", sep="\t", index_col=0)
adata.obs["bac36"] = bac36["norm_score"]


# get Epithelial celltype
Epithelial=adata[adata.obs.category.isin(['Epithelial'])]
sc.pp.filter_cells(Epithelial, min_genes=250)
sc.pp.filter_genes(Epithelial, min_cells=50)
# AnnData object with n_obs × n_vars = 142104 × 21427
# Epithelial.write("/share2/pub/zhenggw/zhenggw/project1/Epithelial.h5ad")

Epithelial=sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/Epithelial.h5ad")
Adult_Epithelial=Epithelial[Epithelial.obs.Age_group.isin(['Adult'])]
#pp
sc.pp.normalize_total(Adult_Epithelial, target_sum=1e4)
sc.pp.log1p(Adult_Epithelial)
sc.pp.highly_variable_genes(Adult_Epithelial, min_mean=0.0125, max_mean=3, min_disp=0.5)
Adult_Epithelial = Adult_Epithelial[:, Adult_Epithelial.var.highly_variable]
sc.pp.regress_out(Adult_Epithelial, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Adult_Epithelial, max_value=10)
sc.tl.pca(Adult_Epithelial, svd_solver='arpack')
sc.pp.neighbors(Adult_Epithelial)

sc.tl.umap(Adult_Epithelial)

#bbknn batch_key = Region code
sc.external.pp.bbknn(Adult_Epithelial, batch_key="Region code")
sc.tl.umap(Adult_Epithelial)
sc.pl.umap(Adult_Epithelial,color='Integrated_05')

#drop the low prop celltype
subEpithelial=Epithelial[Epithelial.obs.Integrated_05.isin(['TA','Enterocyte','Colonocyte','Proximal progenitor','Stem cells','Goblet cell','BEST4+ epithelial','Distal progenitor','Paneth','BEST2+ Goblet cell','Tuft','EECs','Microfold cell','EC cells (TAC1+)'])]

sc.pp.normalize_total(subEpithelial, target_sum=1e4)
sc.pp.log1p(subEpithelial)
sc.pp.highly_variable_genes(subEpithelial, min_mean=0.0125, max_mean=3, min_disp=0.5)
subEpithelial = subEpithelial[:, subEpithelial.var.highly_variable]
sc.pp.regress_out(subEpithelial, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(subEpithelial, max_value=10)
sc.tl.pca(subEpithelial, svd_solver='arpack')
sc.pp.neighbors(subEpithelial)
sc.tl.umap(subEpithelial)

sc.external.pp.bbknn(subEpithelial, batch_key="Region code")
sc.tl.umap(subEpithelial)
sc.pl.umap(subEpithelial,color='Integrated_05')

savefig("Adult_Epithelial_partEpi.pdf")
subEpithelial.write("/share2/pub/zhenggw/zhenggw/project1/subEpithelial.h5ad")



sc.pl.umap(subEpithelial,color="bac27",vmin=0,vmax=5,color_map='RdBu_r')








# get sub_Epithelial celltype
subEpithelial=Epithelial[Epithelial.obs.Integrated_05.isin(['TA','Enterocyte','Colonocyte','Proximal progenitor','Stem cells','Goblet cell','BEST4+ epithelial','Distal progenitor','Paneth','BEST2+ Goblet cell','Tuft','EECs','Microfold cell','EC cells (TAC1+)'])]
sc.pp.filter_cells(subEpithelial, min_genes=250)
sc.pp.filter_genes(subEpithelial, min_cells=50)
# AnnData object with n_obs × n_vars = 140501 × 21296

# subEpithelial.obs.Integrated_05.value_counts()
# TA                     50024
# Enterocyte             37460
# Colonocyte             14133
# Proximal progenitor     9655
# Stem cells              8077
# Goblet cell             5330
# BEST4+ epithelial       4555
# Distal progenitor       3077
# Paneth                  3064
# BEST2+ Goblet cell      3025
# Tuft                     737
# EECs                     485
# Microfold cell           457
# EC cells (TAC1+)         422

#preprocess subEpithelial data
sc.pp.normalize_total(subEpithelial, target_sum=1e4)
sc.pp.log1p(subEpithelial)
sc.pp.highly_variable_genes(subEpithelial, min_mean=0.0125, max_mean=3, min_disp=0.5)
subEpithelial = subEpithelial[:, subEpithelial.var.highly_variable]
sc.pp.regress_out(subEpithelial, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(subEpithelial, max_value=10)
sc.tl.pca(subEpithelial, svd_solver='arpack')
sc.pp.neighbors(subEpithelial)
sc.tl.umap(subEpithelial)
#

