# STOMIX_ST_Analysis_Spot_10x
## Introduction


[STOMIX_ST_Analysis_Spot_10x](the link) pipeline embeds state-of-the-art ST analysis tools, including `Seurat`[1], `sctransform`[2], `Leiden`[3], `BayesSpace`[4], *etc*. to perform spot resolution normalization, imputation, dimension reduction, clustering, differential expression, marker genes enrichment, cell-cell interaction, *etc*.  

## Input files

We require four files from `10x Space Ranger`(https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview) pipelines as input. <mark style="background-color: lightpink">Please directly upload the following files from Space Ranger, no manual amendment is needed. </mark>

1. `Spot expression file` (**required**): the 10x spot expression file with `h5` format, this file is located at `/10x_space_ranger_out_dir/outs` with filename suffix as `filtered_feature_bc_matrix.h5`. 
[Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_filtered_feature_bc_matrix.h5)
The detail description of this file is available at https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/advanced/h5_matrices.

2. `Image file` (**required**): the spatial image, could be high or low resolution with filename suffix as `tissue_lowres_image.png` or `tissue_highres_image.png`, respectively. The image file is located at `/10x_space_ranger_out_dir/outs/spatial`.
![avatar](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_tissue_lowres_image.png)

3. `Spatial coordinates file` (**required**): the spatial coordinate file is located at `/10x_space_ranger_out_dir/outs` with filename suffix as `tissue_positions_list.csv` or `tissue_positions.csv`. 
[Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_tissue_positions_list.csv)
The file columns correspond to the following fields:

    + `barcode`: The sequence of the barcode associated to the spot.
    + `in_tissue`: Binary, indicating if the spot falls inside (1) or outside (0) of tissue.
    + `array_row`: the row id of spot.
    + `array_col`: the column id of spot.
    + `pxl_row_in_fullres`: The row pixel coordinate of the center of the spot in the full resolution image.
    + `pxl_col_in_fullres`: The column pixel coordinate of the center of the spot in the full resolution image.

4. `Spatial scalefactor file` (**required**): the spatial scalefactor json file located in `/10x_space_ranger_out_dir/outs/spatial` with filename suffix as `scalefactors_json.json`.
    [Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_scalefactors_json.json)
    ```
    $ cd /10x_space_ranger_out_dir/outs/spatial
    $ cat scalefactors_json.json
 
    {"tissue_hires_scalef": 0.17011142,
    "tissue_lowres_scalef": 0.051033426",
    "fiducial_diameter_fullres": 144.4773339857,
    "spot_diameter_fullres": 89.43834961021503}
    ```
    This file contains the following fields:
    + `tissue_hires_scalef`: A scaling factor that converts pixel positions in the original, full-resolution image to pixel positions in `tissue_hires_image.png`.
    + `tissue_lowres_scalef`: A scaling factor that converts pixel positions in the original, full-resolution image to pixel positions in `tissue_lowres_image.png`.
    + `fiducial_diameter_fullres`: The number of pixels that span the diameter of a fiducial spot in the original, full-resolution image.
    + `spot_diameter_fullres`: The number of pixels that span the diameter of a theoretical 65µm spot in the original, full-resolution image.


## Input parameters

We require two parameters to fill in.

1. `image_res`: the image resolution you have uploaded. `lowres` for `tisse_lowres_image.png` and `hires` for `tissue_hires_image.png`.

2. `prefix`: the prefix of output files, usually the dataset name. e.g. seurat_v1_hln.


## Output files

The output files are:
- `*_normalized_counts.csv`: the expression counts normalized by `sctransform`. The columns and rows stand for spot barcodes and genes, respectively.
    [Download demo file (head 100 lines).](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_normalized_counts_head100.csv)
    ```
    gene,AAACCGTTCGTCCAGG-1,AAACCTAAGCAGCCGG-1,AAACGAGACGGTTGAT-1
    SOX17, -0.523159629378908,0.921678719367817,-0.80304468337862
    ADHFE1,-0.608935970871144,-1.09977303620335,0.576907392952102
    SULF1,-0.312663321799794,-0.368712005297498,-0.469251763164567
    ```


- `*_spot_meta.csv`: the spot meta information including scaled spatial coordinates, `PCA` and `UMAP` embedding results, `Seurat Leiden` clustering results, *etc*. 
[Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_spot_meta.csv)

    The headers are:
    + `spot_id`: the spot barcode id.
    + `e_UMAP_1`,`e_UMAP_2`: the first two dimensions of UMAP embeddings.
    + `e_PC_1`,`e_PC_2`,`e_PC_3`: the first three dimensions of PCA embeddings.
    + `c_seurat_clusters`: the `Seurat Leiden` clustering results.
    + `c_in_tissue`: binary, indicating if the spot falls inside (1) or outside (0) of tissue.
    + `n_array_row`: the row id of spot.
    + `n_array_col`: the column id of spot.
    + `n_pxl_row_in_fullres`: the row pixel coordinate of the center of the spot in the full resolution image.
    + `n_pxl_col_in_fullres`: the column pixel coordinate of the center of the spot in the full resolution image.
    + `n_pxl_row_in_lowres`: the row pixel coordinate of the center of the spot in the low resolution image.
    + `n_pxl_col_in_lowres`: the column pixel coordinate of the center of the spot in the low resolution image.
    + `n_pxl_row_in_hires`: the row pixel coordinate of the center of the spot in the high resolution image.
    + `n_pxl_col_in_hires`: the column pixel coordinate of the center of the spot in the high resolution image.

    Please note that prefix `n` denotes numeric value, prefix `c` denotes category meta information like group or cluster. Prefix `e` refers to embedding/dimension reduction methods.

- `*_marker_genes.csv`: the list of marker genes. 
    [Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_marker_genes.csv)
    ```
    CYP2A5
    SLC22A6
    CYP2E1
    PCK1
    SPINK1
    ```

- `*_de_markers.csv`: the differential expressed markers among clusters detected by `Seurat FindAllMarkers` (https://satijalab.org/seurat/reference/findallmarkers). 
    [Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_de_markers.csv)
    ```
    gene,p_val,avg_log2FC,pct.1,pct.2,p_val_adj,c_seurat_clusters
    CYP2A5,2.14852567803993e-123,2.12438233188317,1,0.853,3.67333435174486e-119,0
    SLC22A6,8.67255144425018e-121,1.73364344386156,1,0.805,1.48274612042345e-116,0
    CYP2E1,1.84478367314248e-118,1.82951971472933,1,0.939,3.15402664597169e-114,0
    PCK1,3.00294594308557e-118,1.49509186919739,1,0.996,5.1341366788934e-114,0
    SPINK1,1.95750838346582e-117,1.0171935406928,1,0.999,3.34675208321152e-113,0
    ```
    The file headers are:
    + `gene`: the Hugo gene name.
    + `p_val`: the p value of Wilcoxon Rank Sum test.
    + `avg_log2FC`: log fold-chage of the average expression among clusters.
    + `pct.1`: the percentage of spots where the feature is detected inside the specified cluster.
    + `pct.2`: the percentage of spots where the feature is detected outside the specified cluster.
    + `p_val_adj`: adjusted p-value, based on bonferroni correction using all features in the dataset.

- `*_marker_gene_exp.csv`: the expression profiles of marker genes. 
    [Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_marker_gene_exp.csv)
    ```
    IGFBP2,PAM,TFCP2L1,PIGR,CTSE
    -0.200407506774145,1.84930980776098,0.984699256618785,0.0141368340887279,-0.231425061669138
    -0.24981065085656,-0.272777391877713,0.447035163685717,4.55161252519223,-0.281968753930634
    -0.338726119548364,-0.292249069452559,-0.571216723814115,-1.31179674036995,-0.372865114113906
    ```
- `*_marker_gene_stats.csv`: the quantile stats of marker genes from `Seurat Leiden` cluster. 
    [Download demo file.](https://raw.githubusercontent.com/deepomicslab/STOMIX/main/demo_data/10xDemoMK_mouse-kidney-section-coronal-1-standard-1-1-0_section1_slice1_marker_gene_stat.csv)
    ```
    cluster,gene,mean,min,Q1,median,Q3,max
    2,IGFBP2,0.0242452513738838,-0.267623361545279,-0.189775924340351,-0.1649898401019,-0.132699041924791,6.18999278998317
    2,PAM,0.259606538357579,-1.46668363707726,-0.746182046947719,0.0211419388766112,1.14945269587931,4.5906233753646
    2,TFCP2L1,0.14705664357956,-1.82812381613413,-0.590682381481341,-0.109587125367983,0.78553241438755,3.14880717654735
    ```
    The file headers are:
    + `cluster`: the `Seurat Leiden` cluster id.
    + `gene`: the differential expressed gene name.
    + `mean`: the mean expression of gene in cluster.
    + `min`: the minimal expression of gene in cluster.
    + `Q1`: the 25% quantile expression of gene in cluster.
    + `median`: the median expression of gene in cluster.    
    + `Q3`: the 75% quantile expression of gene in cluster.
    + `max`: the maximal expression of gene in cluster. 

## Output Visualizations


## References

[1] Satija, R., Farrell, J.A., Gennert, D., Schier, A.F. and Regev, A., 2015. Spatial reconstruction of single-cell gene expression data. *Nature biotechnology*, 33(5), pp.495-502.

[2] Hafemeister, C. and Satija, R., 2019. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome biology*, 20(1), p.296.

[3] Traag, V.A., Waltman, L. and Van Eck, N.J., 2019. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific reports*, 9(1), p.5233.

[4] Zhao, E., Stone, M.R., Ren, X., Guenthoer, J., Smythe, K.S., Pulliam, T., Williams, S.R., Uytingco, C.R., Taylor, S.E., Nghiem, P. and Bielas, J.H., 2021. Spatial transcriptomics at subspot resolution with BayesSpace. *Nature Biotechnology*, 39(11), pp.1375-1384.




