```R
library(Seurat)
library(SeuratDisk)

Convert("20210204_S_SUBS4_4x_annotated_donors_UMI_counts_HCA.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc <- LoadH5Seurat("20210204_S_SUBS4_4x_annotated_donors_UMI_counts_HCA.h5seurat")
pbmc

test <- CreateSeuratObject(pbmc@assays$RNA@counts,project="test")
pbmc@meta.data$nCount_RNA <- test@meta.data$nCount_RNA
pbmc@meta.data$nFeature_RNA <- test@meta.data$nFeature_RNA
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
```

```R
library(cowplot)
library(harmony)
library(ggplot2)

pbmc <- RunHarmony(pbmc,group.by.vars="donor",plot_convergence = TRUE)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20) %>% FindClusters()
saveRDS(pbmc, file = "./pbmc.rds")
save(pbmc,file = "./pbmc.RData")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10,file="top10markers.csv")

pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC) -> top2
VlnPlot(pbmc,features=top2$gene,stack=T)

markers <- top2$gene[!duplicated(top2$gene)]
VlnPlot(pbmc,features=markers,stack=T)+NoLegend()
DotPlot(pbmc,features=markers)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

0	CEBPB		HSC/MPP1
1	KLF6		HSC/MPP2
2	VIM			HSC/MPP3
3	SERPINB1		MEP1
4	CA1				EryP1
5	HIST1H1C	HSC/MPP4
6	LTB				LP1
7	C1QTNF4		HSC/MPP5
8	KLF6		HSC/MPP6
9	TYMS			MP1
10	RANBP1			MEP2
11	VIM			HSC/MPP7
12	TYMS			MP3
13	IGLL1			MP4
14	MT-ATP6			MP5
15	MPO				MP6
16	CNRIP1			MEP3
17	IRF8			DC
18	DNAJB1		HSC/MPP8
19	AVP			HSC/MPP9
20	CENPF	  HSPC cycle1
21	MS4A3			BCMP
22	ITGA2B	  		MkP
23	FAM211A		HSC/MPP10
24	MKI67	  HSPC cycle2
25	MT-CO3		   EryP2
26	MALAT1		    LP2


new.cluster.ids <- c("HSC/MPP1", "HSC/MPP2", "HSC/MPP3", "MEP1", "EryP1", "HSC/MPP4","LP1","HSC/MPP5","HSC/MPP6","MP1","MEP2","HSC/MPP7","MP3","MP4",
    "MP5", "MP6", "MEP3", "DC", "HSC/MPP8", "HSC/MPP9", "HSPC cycle1", "BCMP", "MkP", "MEP4", "HSPC cycle2", "EryP2", "LP2")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc <- StashIdent(pbmc, save.name = 'new.cluster.ids')

save(pbmc,file="pbmc.RData")


pbmc.loom <- as.loom(pbmc,file="pbmc.loom",verbose = FALSE)
pbmc.loom
pbmc.loom$close_all()


library(SeuratDisk)

SaveH5Seurat(pbmc,filename="3site.h5Seurat")
Convert("3site.h5Seurat",dest="h5ad")

#/share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/
```

```python
import scanpy as sc
import pandas as pd
import os
#import loompy
import matplotlib.pyplot as plt

raw_adata = sc.read_h5ad("3site.h5ad")
raw_adata.obs['leiden'] = raw_adata.obs['new.cluster.ids'].astype("category")

#add celltype
new_cluster_names = ['HSC/MPP1', 'HSC/MPP2', 'HSC/MPP3', 'MEP1', 'EryP1', 'HSC/MPP4','LP1','HSC/MPP5','HSC/MPP6','MP1','MEP2','HSC/MPP7','MP3','MP4',
    'MP5', 'MP6', 'MEP3', 'DC', 'HSC/MPP8', 'HSC/MPP9', 'HSPC cycle1', 'BCMP', 'MkP', 'MEP4', 'HSPC cycle2', 'EryP2', 'LP2']
raw_adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(raw_adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')


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

```shell
#!/bin/sh
#PBS -N AdultHemOrgans
#PBS -q workq
#PBS -l mem=100G
#PBS -o /share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/AdultHemOrgans.out
#PBS -e /share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/AdultHemOrgans.err

/share2/pub/zhenggw/zhenggw/anaconda3/envs/python3/bin/python /share2/pub/zhenggw/zhenggw/scDRS/compute_score.py \
    --h5ad_file /share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/expr.h5ad \
    --h5ad_species human \
    --gs_file /share2/pub/lijj/lijj/TheDutchMicrobiomeProject/GenusLevel/scDRS/Bac_munge.gs \
    --gs_species human \
    --cov_file /share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/cov.tsv \
    --flag_filter True \
    --flag_raw_count True \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score True \
    --out_folder /share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/data/ > /share2/pub/zhenggw/zhenggw/project1/AdultHemOrgans/scDRS.log 2>&1
```

