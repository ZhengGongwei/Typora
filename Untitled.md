

```python
import scdrs
from scipy import stats
import pandas as pd
import scanpy as sc
from anndata import AnnData
sc.set_figure_params(dpi=125)
import numpy as np
import statsmodels.api as sm

import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings("ignore")

os.chdir('/share2/pub/zhenggw/zhenggw/project1/test/')
adata = sc.read_h5ad("data/expr.h5ad")
df_gs = pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/test/data/geneset.gs", sep="\t", index_col=0)

dict_score = {
    trait: pd.read_csv(f"data/{trait}.full_score.gz", sep="\t", index_col=0)
    for trait in df_gs.index
}


sc.set_figure_params(figsize=[1, 1], dpi=150)
plt.figure(dpi=300,figsize=(24,8))
plt.xticks(fontsize=10)

#dict() creates a dictionary
dict_stats = dict()
for trait in dict_score:
    dict_stats[trait] = scdrs.method.group_stats(
        dict_score[trait], adata, group_col="level1class"
    )
 
with open('test.csv', 'w') as f:
    for key in dict_stats.keys():
        f.write("%s, %s\n" % (key, dict_stats[key]))


with open('maintest.csv', 'w', newline='') as csvfile:
    header_key = ['trait', 'data']
    new_val = csv.DictWriter(csvfile, fieldnames=header_key)
    new_val.writeheader()
    for new_k in dict_stats.keys():
        new_val.writerow({'trait': new_k, 'data': dict_stats[new_k]})  

 dict_stats.to_csv("data/dict_stats.tsv", sep="\t")

np.save('data/dict_stats.npy', dict_stats) # 注意带上后缀名

load_dict = np.load('data/dict_stats.npy',allow_pickle=True).item()

#visualize the  statistics simulateously using the built-in scdrs.util.plot_group_stats

dict_celltype_display_name = {
    "pyramidal_CA1": "Pyramidal CA1",
    "oligodendrocytes": "Oligodendrocyte",
    "pyramidal_SS": "Pyramidal SS",
    "interneurons": "Interneuron",
    "endothelial-mural": "Endothelial",
    "astrocytes_ependymal": "Astrocyte",
    "microglia": "Microglia",
}

scdrs.util.plot_group_stats(
    {
        trait: df_stats.rename(index=dict_celltype_display_name)
        for trait, df_stats in dict_stats.items()
    }
)

plt.savefig("filename.png")
plt.show()
```