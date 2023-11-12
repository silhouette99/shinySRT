## H5ad process

h5ad is a common way of storing data in python, which can be called directly without memory, and scanpy is a very mature spatial transcriptome data processing tool in python.
The following is still using 10x Visim public data [Mouse Brain Serial Section 2 (Sagittal-Anterior)](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard), including the data h5ad transformation and `ShinySRT` example code:


``` python
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
import os
os.chdir('./')
results_file = 'Anterior.h5ad' 

## data storied into the directory
adata = sc.read_visium('anterior/')
adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.match("^MT-|^mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata.var["RP"] = adata.var_names.str.match("^RPS|^RPL|^Rps|^Rpl")
sc.pp.calculate_qc_metrics(adata, qc_vars=["RP"], inplace=True)

sc.pp.filter_cells(adata, min_counts=5000) #filtered out 80 cells that have less than 5000 counts
sc.pp.filter_cells(adata, max_counts=50000) #filtered out 39 cells that have more than 50000 counts
adata = adata[adata.obs["pct_counts_mt"] < 25]
print(f"#cells after MT filter: {adata.n_obs}") #cells after MT filter: 2502
sc.pp.filter_genes(adata, min_cells=10)


sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)


sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")

adata.write(results_file)

```

``` r
makespashiny(dat = 'Anterior.h5ad',title = 'spatial experiment',gex.assay = 'counts')
```