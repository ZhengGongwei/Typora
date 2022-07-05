```R
library(Seurat)
covid <- readRDS("../COVID_combined/covid_combined_annotation.rds")
head(covid@meta.data)
                                 orig.ident nCount_RNA nFeature_RNA percent.mt   RNA_snn_res.0.5   seurat_clusters     annotation
    AAACGGGCATCCCACT-1-1:1_1:1-1  Covid_46W   6986.000         1167   1.535756       6                  6   	       Naive CD8+T
    AAATGCCCATATACCG-1-1:1_1:1-1  Covid_46W   4208.530          609   2.337488       5                  5  			   Memory CD8+T
covid
    An object of class Seurat
    18092 features across 1116364 samples within 1 assay
    Active assay: RNA (18092 features, 2000 variable features)
    3 dimensional reductions calculated: pca, harmony, umap

library(SeuratDisk)
SaveH5Seurat(covid,filename="/share2/pub/zhenggw/zhenggw/covid100W_scRDS/covid.h5Seurat")
Convert("/share2/pub/zhenggw/zhenggw/covid100W_scRDS/covid.h5Seurat",dest="h5ad")
```

```python
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

raw_adata = sc.read_h5ad("covid.h5ad")
raw_adata.obs['annotation'] = raw_adata.obs['annotation'].astype("category")

adata = raw_adata.raw.to_adata()
adata.obsp = raw_adata.obsp
adata.X = adata.X

# compile covariates
df_cov = pd.DataFrame(index=adata.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = adata.obs["nFeature_RNA"]

if not os.path.exists("data/"):
	os.makedirs("data/")

adata.write_h5ad("data/expr.h5ad")
df_cov.to_csv("data/cov.tsv", sep="\t")
```

```R
library(gprofiler2)
library(readr)
data <- read.table("magma.tsv",header=T)
a <- gconvert(data$GENE,organism="hsapiens",numeric_ns="ENTREZGENE_ACC",target="ENSG",mthreshold = Inf, filter_na = TRUE)
aa <- merge(data,a,by.x="GENE",by.y="input")
gene <- aa[,c(6,2)]
colnames(gene) <- c("GENE","COVID")
write_tsv(gene,file="covidgene.tsv")
```

```shell
scdrs munge-gs \
    --out-file covid.gs \
    --zscore-file covidgene.tsv \
    --weight zscore  \
    --n-max 1000
 
 
 sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
```

