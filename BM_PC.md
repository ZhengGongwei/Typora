```R
library(Seurat)
load("BMPC_merge.Rdata")
BMPC <- sce.big
rm(sce.big)
BMPC@meta.data$orig.ident <- "BMPC"

BMPC[["percent.mt"]] <- PercentageFeatureSet(BMPC, pattern = "^MT-")
metadata <- read.table("GSE117156_metadata.txt.gz",header=T)
BMPC@meta.data$Amp_batch_ID <- metadata$Amp_batch_ID

id <- metadata$Experiment_ID
for(i in 1:length(id)){
    id[i] <- substring(metadata$Experiment_ID[i],9)
}

BMPC@meta.data$Experiment_ID <- id

save(BMPC,file="BMPC.RData")

# Visualize QC metrics as a violin plot
Idents(BMPC) <- "orig.ident"
VlnPlot(BMPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

BMPC <- subset(BMPC, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 50)

BMPC <- NormalizeData(BMPC)

BMPC <- FindVariableFeatures(BMPC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BMPC)
BMPC <- ScaleData(BMPC, features = all.genes)
```



