```python
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

adata = sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True)
#save 1counts.png
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

#Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
#save 2pca.png

#Computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color='batch')
```

```python
#use bbknn
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import bbknn

adata = sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch', 'category'], palette=sc.pl.palettes.vega_20_scanpy)
sc.pl.umap(adata, color='category', palette=sc.pl.palettes.vega_20_scanpy)


#bbknn
sc.external.pp.bbknn(adata, batch_key='batch')
sc.tl.umap(adata)
sc.pl.umap(adata, color='category')
adata.write_h5ad("expr_RB.h5ad")

#subset
adata.obs.Diagnosis.value_counts()
fetal                      247344
Healthy adult              124367
Pediatric healthy           29594
Pediatric Crohn Disease     27164
Name: Diagnosis, dtype: int64

adata.obs.Age_group.value_counts()
First trim         149224
Adult              113203
Second trim         89946
Pediatric           29594
Pediatric_IBD       27164
Adult_MLN           11164
Second trim_MLN      8174
Name: Age_group, dtype: int64

        
Name: Age_group, Length: 428469, dtype: category
Categories (7, object): ['Adult', 'Adult_MLN', 'First trim', 'Pediatric', 'Pediatric_IBD',
                         'Second trim', 'Second trim_MLN']

Adult=adata[adata.obs.Age_group.isin(['Adult'])]
#Adult.write_h5ad("Adult.h5ad")
Adult_MLN=adata[adata.obs.Age_group.isin(['Adult_MLN'])]
#Adult_MLN.write_h5ad("Adult_MLN.h5ad")
First_trim=adata[adata.obs.Age_group.isin(['First trim'])]
#First_trim.write_h5ad("First_trim.h5ad")
Pediatric=adata[adata.obs.Age_group.isin(['Pediatric'])]
#Pediatric.write_h5ad("Pediatric.h5ad")
Pediatric_IBD=adata[adata.obs.Age_group.isin(['Pediatric_IBD'])]
#Pediatric_IBD.write_h5ad("Pediatric_IBD.h5ad")
Second_trim=adata[adata.obs.Age_group.isin(['Second trim'])]
#Second_trim.write_h5ad("Second_trim.h5ad")
Second_trim_MLN=adata[adata.obs.Age_group.isin(['Second trim_MLN'])]
#Second_trim_MLN.write_h5ad("Second_trim_MLN.h5ad")


import matplotlib.pyplot as plt
 
labels = ['Epithelial', 'Plasma cells',"T cells","Mesenchymal","Endothelial","Myeloid","B cells","Neuronal"]

Adult_number = Adult.obs.category.value_counts()
Adult_MLN_number = Adult_MLN.obs.category.value_counts()
First_trim_number = First_trim.obs.category.value_counts()
Pediatric_number = Pediatric.obs.category.value_counts()
Pediatric_IBD_number = Pediatric_IBD.obs.category.value_counts()
Second_trim_number = Second_trim.obs.category.value_counts()
Second_trim_MLN_number = Second_trim_MLN.obs.category.value_counts()
frames = [Adult_number, Adult_MLN_number, First_trim_number, Pediatric_number, Pediatric_IBD_number, Second_trim_number, Second_trim_MLN_number]
result = pd.concat(frames,axis=1,keys=["Adult", "Adult_MLN", "First_trim", "Pediatric", "Pediatric_IBD", "Second_trim", "Second_trim_MLN"])
result.fillna (0,inplace=True)
cell_prop=result.div(result.sum(axis=0), axis=1)
```

```R
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(RColorBrewer)
library(reshape2)
result <- read.table("celltype_prop.txt",header = T,sep="\t")
result

result2 <- melt(result,id.vars='celltype')
palette = brewer.pal(9,'RdYlBu')
ggplot(data=result2,aes(variable,value,fill=celltype))+
  geom_bar(stat="identity", position="fill", width=0.8,size=0.25)+
  scale_fill_manual(values=palette)+
  xlab("Samples") + ylab("Expression values")+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    legend.title=element_text(size=14,face="plain",color="black"),
    legend.position = "right"
  	)
```

```python
adata=sc.read_h5ad("expr_RB.h5ad")
Adult=adata[adata.obs.Age_group.isin(['Adult'])]
Adult_MLN=adata[adata.obs.Age_group.isin(['Adult_MLN'])]
First_trim=adata[adata.obs.Age_group.isin(['First trim'])]
Pediatric=adata[adata.obs.Age_group.isin(['Pediatric'])]
Pediatric_IBD=adata[adata.obs.Age_group.isin(['Pediatric_IBD'])]
Second_trim=adata[adata.obs.Age_group.isin(['Second trim'])]
Second_trim_MLN=adata[adata.obs.Age_group.isin(['Second trim_MLN'])]

sc.pl.umap(Adult,color='category')
sc.pl.umap(Adult_MLN,color='category')
sc.pl.umap(First_trim,color='category')
sc.pl.umap(Pediatric,color='category')
sc.pl.umap(Pediatric_IBD,color='category')
sc.pl.umap(Second_trim,color='category')
sc.pl.umap(Second_trim_MLN,color='category')

df_cov = pd.DataFrame(index=Adult.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = Adult.obs["n_genes"]
df_cov.to_csv("Adult_cov.tsv", sep="\t")

df_cov = pd.DataFrame(index=Adult_MLN.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = Adult_MLN.obs["n_genes"]
df_cov.to_csv("Adult_MLN_cov.tsv", sep="\t")

df_cov = pd.DataFrame(index=First_trim.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = First_trim.obs["n_genes"]
df_cov.to_csv("First_trim_cov.tsv", sep="\t")

df_cov = pd.DataFrame(index=Pediatric.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = Pediatric.obs["n_genes"]
df_cov.to_csv("Pediatric_cov.tsv", sep="\t")

df_cov = pd.DataFrame(index=Pediatric_IBD.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = Pediatric_IBD.obs["n_genes"]
df_cov.to_csv("Pediatric_IBD_cov.tsv", sep="\t")

df_cov = pd.DataFrame(index=Second_trim.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = Second_trim.obs["n_genes"]
df_cov.to_csv("Second_trim_cov.tsv", sep="\t")

df_cov = pd.DataFrame(index=Second_trim_MLN.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = Second_trim_MLN.obs["n_genes"]
df_cov.to_csv("Second_trim_MLN_cov.tsv", sep="\t")
```

