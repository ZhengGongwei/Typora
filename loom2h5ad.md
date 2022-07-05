```python
import scanpy as sc
import pandas as pd
import loompy
import os
#os.getwd()
#os.listdir()

adata_ref = sc.read_loom("./6d046877-2048-471e-bd3b-ea198a493022.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
adata_ref.var_names_make_unique()
adata1 = adata_ref

adata_ref = sc.read_loom("./progenitor-hierarchy-mouse-adipose-10XV2.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
```

