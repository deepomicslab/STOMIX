# TIMEDB cell Correlation
# Introduction
With TIME cell fraction as input, we conducted the Pearson correlation and associated significant test with R package “corrplot”. P value less than 0.05 is regarded as statistically significant.


### Citation
 Taiyun Wei, Viliam Simko, Michael Levy, Yihui Xie, Yan Jin, and Jeff Zemla. Package ‘cor-rplot’ Statistician, 56(316):e24, 2017.
  
## Module Structure


## Input file
`Cell fraction data`(**required**):the TIME cell fraction result(csv).
*  The first two columns in the first row should be "sample_name" and "method", followed by the TIME cell type
*  The first column should be the the sample name.
*  The second column should be the method name to get the TIME cell fraction result.
![avatar](./plots/input_cell_fraction.jpg)

`Cluster result data`(**optional**):the cluster result file(csv)
  - First row should be a header with a 'sample_name'column lable followed by clinical features.
  - The header line keyword of classfication feature follows c_[feature] format, such as 'c_gender'.
  - The header line keyword of continuous feature follows n_[feature] format, such as 'c_age'.
  - The column related to survival days should be named "os" or "pfs"(**required**), details are:
    - os: available for cancer or diease donor, overall survival days, the length of time from either the date of diagnosis or the start of treatment, that a patient still alive.
    - pfs: available for cancer or diease donor, progression-free survival days, the length of time during and after the treatment, that a patient lives with the disease but it does not get worse.
    - os_status: available for cancer donor, the os outcome, binary, value of 1 for death, 0 for alive.
    - pfs_status: available for cancer donor, the pfs outcome, binary, value of 1 for pregression or recurrence, 0 for otherwise.
  - More details could be seen in our demo file.
![avatar](./plots/input_clinical_data.jpg)

## Results

### Correlation result file(csv)
The correlation result file: We classify sample according to "c_" column of cluster result file, and analyze correlation between fraction values of two different immune cell types. Value column is correlation value, pvalue is significant value, metric column is correlation method, method is cell fraction analysis method, feature and group indicate current classification of sample.

![avatar](./plots/correlation.jpg)
### Visualization





## Download demo data
- Cell fraction result file(./data/TCGA_ACC_quanTIseq.csv)
- Cluster result file(./data/quanTIseq_cluster_Result.csv)
- Correlation result file(./data/demo_Correlation.csv)