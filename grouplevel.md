```python
import scdrs
from scipy import stats
import pandas as pd
import scanpy as sc
from anndata import AnnData
#sc.set_figure_params(dpi=125)
import numpy as np
import statsmodels.api as sm

import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings("ignore")

os.chdir('/share2/pub/zhenggw/zhenggw/project1/intest2/Pediatric/')
Pediatric = sc.read_h5ad("Pediatric.h5ad")
df_gs = pd.read_csv("/share2/pub/zhenggw/zhenggw/project1/intest/umap/Bac_munge.gs", sep="\t", index_col=0)

dict_score = {
    trait: pd.read_csv(f"output_RB/{trait}.full_score.gz", sep="\t", index_col=0)
    for trait in df_gs.index
}


dict_stats = dict()
for trait in dict_score:
    dict_stats[trait] = scdrs.method.group_stats(
        dict_score[trait], Pediatric, group_col="category"
    )
    
np.save('data/dict_stats.npy', dict_stats) # 注意带上后缀名
    
with open('intest.csv', 'w') as f:
    for key in dict_stats.keys():
        f.write("%s, %s\n" % (key, dict_stats[key]))

#visualize the  statistics simulateously using the built-in scdrs.util.plot_group_stats
# dict_celltype_display_name = {
#     'B cells': "B",
#     'T cells': "T",
#     'Epithelial': "Epithelial",
#     'Mesenchymal': "Epithelial",
#     'Plasma cells': "Plasma",
#     'Myeloid': "Myeloid",
#     'Endothelial': "Endothelial",
#     'Neuronal': "Neuronal",
#     'Red blood cells': "Red blood",
# }

# scdrs.util.plot_group_stats(
#     {
#         trait: df_stats.rename(index=dict_celltype_display_name)
#         for trait, df_stats in dict_stats.items()
#     }
# )

# plt.savefig("IntestGrouplevel.png")
```

