# STOMIX_ST_10x_Spot_Analysis
## Introduction


[STOMIX_ST_10x_Spot_Analysis](the link) pipeline embeds state-of-the-art ST analysis tools, including `Seurat`, `sctransform`, `Leiden`, *etc*. to perform spot resolution normalization, imputation, dimension reduction, clustering, differential expression, marker genes enrichment, cell-cell interaction, *etc*.  


## Input

It requires four input files from `10x Space Ranger` and two parameters:

```
Input files:
 --in_h5_fn=The path of 10x ST h5 file.                 e.g. filtered_feature_bc_matrix.h5
 --in_image_fn=The path of spatial image.               e.g. tissue_lowres_image.png
 --in_coord_fn=The path of 10x spatial coordinates.     e.g. tissue_positions_list.csv
 --in_scale_fn=The path of scalefactor json.            e.g. scalefactors_json.json

Parameters:
 --image_res=Can only be "lowres" or "hires".           e.g. lowres
 --prefix=The prefix of output file.                    e.g. seurat_v1_hln

```

Download demo files:
- [in_h5_fn]()
- [in_image_fn]()
- [in_coord_fn]()
- [in_scale_fn]()

## Output
The output files are:
- `*_normalized_counts.csv`: the expression counts normalized by `sctransform`. [Download demo.]()
- `*_spot_meta.csv`: the spot meta information including scaled spatial coordinates, `PCA` and `UMAP` embedding results, `Seurat Leiden` clustering results, *etc*. [Download demo.]()
- `*_marker_genes.csv`: the list of marker genes. [Download demo.]()
- `*_de_markers.csv`: the differential expressed markers among detected clusters. [Download demo.]()
- `*_marker_gene_exp.csv`: the expression profiles of marker genes. [Download demo.]()
- `*_marker_gene_stats.csv`: the quantile stats of marker genes from `Seurat Leiden` cluster.

  
## Visualization



