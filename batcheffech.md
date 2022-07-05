```python
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

adata = sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
MLN=adata[adata.obs.Age_group.isin(['Adult_MLN'])]

sc.pp.normalize_total(MLN, target_sum=1e4)
sc.pp.log1p(MLN)
sc.pp.filter_genes(MLN, min_cells=3)

sc.pp.regress_out(MLN, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(MLN, max_value=10)
sc.tl.pca(MLN, svd_solver='arpack')
sc.pp.neighbors(MLN, n_neighbors=10, n_pcs=40)
sc.tl.umap(MLN)
sc.pl.umap(MLN,color='category')

sc.external.pp.bbknn(MLN, batch_key='batch')
sc.tl.umap(MLN)
sc.pl.umap(MLN,color='category')
MLN.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/MLN_BBKNN/MLN.h5ad")

df_cov = pd.DataFrame(index=MLN.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = MLN.obs["n_genes"]
df_cov.to_csv("/share2/pub/zhenggw/zhenggw/project1/intest2/MLN_BBKNN/MLN_cov.tsv", sep="\t")
```

```python
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt

adata = sc.read_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest/Full_obj_raw_counts_nosoupx.h5ad")
#1
Adult=adata[adata.obs.Age_group.isin(['Adult'])]
#Adult.write_h5ad("Adult.h5ad")

sc.pp.normalize_total(Adult, target_sum=1e4)
sc.pp.log1p(Adult)
sc.pp.filter_genes(Adult, min_cells=3)

sc.pp.regress_out(Adult, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Adult, max_value=10)
sc.tl.pca(Adult, svd_solver='arpack')
sc.pp.neighbors(Adult, n_neighbors=10, n_pcs=40)
sc.tl.umap(Adult)
sc.pl.umap(Adult,color='category')
plt.savefig("Adult.png")
plt.show()

sc.external.pp.bbknn(Adult, batch_key='batch')
sc.tl.umap(Adult)
sc.pl.umap(Adult,color='category')
plt.savefig("Adult_bbknn.png")
plt.show()

Adult.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/Adult/Adult.h5ad")

#2
Adult_MLN=adata[adata.obs.Age_group.isin(['Adult_MLN'])]
#Adult_MLN.write_h5ad("Adult_MLN.h5ad")

sc.pp.normalize_total(Adult_MLN, target_sum=1e4)
sc.pp.log1p(Adult_MLN)
sc.pp.filter_genes(Adult_MLN, min_cells=3)

sc.pp.regress_out(Adult_MLN, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Adult_MLN, max_value=10)
sc.tl.pca(Adult_MLN, svd_solver='arpack')
sc.pp.neighbors(Adult_MLN, n_neighbors=10, n_pcs=40)
sc.tl.umap(Adult_MLN)
sc.pl.umap(Adult_MLN,color='category')
plt.savefig("Adult_MLN.png")
plt.show()

sc.external.pp.bbknn(Adult_MLN, batch_key='batch')
sc.tl.umap(Adult_MLN)
sc.pl.umap(Adult_MLN,color='category')
plt.savefig("Adult_MLN_bbknn.png")
plt.show()

Adult_MLN.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/Adult_MLN/Adult_MLN.h5ad")

#3
First_trim=adata[adata.obs.Age_group.isin(['First trim'])]
#First_trim.write_h5ad("First_trim.h5ad")
sc.pp.normalize_total(First_trim, target_sum=1e4)
sc.pp.log1p(First_trim)
sc.pp.filter_genes(First_trim, min_cells=3)

sc.pp.regress_out(First_trim, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(First_trim, max_value=10)
sc.tl.pca(First_trim, svd_solver='arpack')
sc.pp.neighbors(First_trim, n_neighbors=10, n_pcs=40)
sc.tl.umap(First_trim)
sc.pl.umap(First_trim,color='category')
plt.savefig("First_trim.png")
plt.show()

sc.external.pp.bbknn(First_trim, batch_key='batch')
sc.tl.umap(First_trim)
sc.pl.umap(First_trim,color='category')
plt.savefig("First_trim_bbknn.png")
plt.show()

First_trim.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/First_trim/First_trim.h5ad")

#4
Pediatric=adata[adata.obs.Age_group.isin(['Pediatric'])]
#Pediatric.write_h5ad("Pediatric.h5ad")

sc.pp.normalize_total(Pediatric, target_sum=1e4)
sc.pp.log1p(Pediatric)
sc.pp.filter_genes(Pediatric, min_cells=3)

sc.pp.regress_out(Pediatric, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Pediatric, max_value=10)
sc.tl.pca(Pediatric, svd_solver='arpack')
sc.pp.neighbors(Pediatric, n_neighbors=10, n_pcs=40)
sc.tl.umap(Pediatric)
sc.pl.umap(Pediatric,color='category')
plt.savefig("Pediatric.png")
plt.show()

sc.external.pp.bbknn(Pediatric, batch_key='batch')
sc.tl.umap(Pediatric)
sc.pl.umap(Pediatric,color='category')
plt.savefig("Pediatric_bbknn.png")
plt.show()

Pediatric.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/Pediatric/Pediatric.h5ad")
#5
Pediatric_IBD=adata[adata.obs.Age_group.isin(['Pediatric_IBD'])]
#Pediatric_IBD.write_h5ad("Pediatric_IBD.h5ad")
sc.pp.normalize_total(Pediatric_IBD, target_sum=1e4)
sc.pp.log1p(Pediatric_IBD)
sc.pp.filter_genes(Pediatric_IBD, min_cells=3)

sc.pp.regress_out(Pediatric_IBD, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Pediatric_IBD, max_value=10)
sc.tl.pca(Pediatric_IBD, svd_solver='arpack')
sc.pp.neighbors(Pediatric_IBD, n_neighbors=10, n_pcs=40)
sc.tl.umap(Pediatric_IBD)
sc.pl.umap(Pediatric_IBD,color='category')
plt.savefig("Pediatric_IBD.png")
plt.show()

sc.external.pp.bbknn(Pediatric_IBD, batch_key='batch')
sc.tl.umap(Pediatric_IBD)
sc.pl.umap(Pediatric_IBD,color='category')
plt.savefig("Pediatric_IBD_bbknn.png")
plt.show()

Pediatric_IBD.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/Pediatric_IBD/Pediatric_IBD.h5ad")

#6
Second_trim=adata[adata.obs.Age_group.isin(['Second trim'])]
#Second_trim.write_h5ad("Second_trim.h5ad")
sc.pp.normalize_total(Second_trim, target_sum=1e4)
sc.pp.log1p(Second_trim)
sc.pp.filter_genes(Second_trim, min_cells=3)

sc.pp.regress_out(Second_trim, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Second_trim, max_value=10)
sc.tl.pca(Second_trim, svd_solver='arpack')
sc.pp.neighbors(Second_trim, n_neighbors=10, n_pcs=40)
sc.tl.umap(Second_trim)
sc.pl.umap(Second_trim,color='category')
plt.savefig("Second_trim.png")
plt.show()

sc.external.pp.bbknn(Second_trim, batch_key='batch')
sc.tl.umap(Second_trim)
sc.pl.umap(Second_trim,color='category')
plt.savefig("Second_trim_bbknn.png")
plt.show()

Second_trim.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/Second_trim/Second_trim.h5ad")

#7
Second_trim_MLN=adata[adata.obs.Age_group.isin(['Second trim_MLN'])]
#Second_trim_MLN.write_h5ad("Second_trim_MLN.h5ad")
sc.pp.normalize_total(Second_trim_MLN, target_sum=1e4)
sc.pp.log1p(Second_trim_MLN)
sc.pp.filter_genes(Second_trim_MLN, min_cells=3)

sc.pp.regress_out(Second_trim_MLN, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(Second_trim_MLN, max_value=10)
sc.tl.pca(Second_trim_MLN, svd_solver='arpack')
sc.pp.neighbors(Second_trim_MLN, n_neighbors=10, n_pcs=40)
sc.tl.umap(Second_trim_MLN)
sc.pl.umap(Second_trim_MLN,color='category')
plt.savefig("Second_trim_MLN.png")
plt.show()

sc.external.pp.bbknn(Second_trim_MLN, batch_key='batch')
sc.tl.umap(Second_trim_MLN)
sc.pl.umap(Second_trim_MLN,color='category')
plt.savefig("Second_trim_MLN_bbknn.png")
plt.show()

Second_trim_MLN.write_h5ad("/share2/pub/zhenggw/zhenggw/project1/intest2/Second_trim_MLN/Second_trim_MLN.h5ad")
```

