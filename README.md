# Scripts for WGCNA usage
To preprocess an expression matrix and perform Weighted Gene Co-expression Network Analysis (WGCNA). The preprocessing steps include removing low expression genes, selecting focal tissue samples, and filtering genes based on their coefficient of variation. The WGCNA script constructs a gene co-expression network, identifies gene modules, and relates them to external traits.

## Data Preprocessing

The preprocessing script `01_datExpr_preprocessing.ipynb` performs the following steps:

1. Remove Low Expression Genes:

Filters out genes based on a TPM cutoff.
Plots the number of genes passing different TPM thresholds.

2. Select Focal Tissue/Organ Samples:

Filters the expression matrix to include only samples from a specified tissue/organ.

3. Remove Housekeeping Genes:

Filters out housekeeping genes based on their coefficient of variation (CV).
Plots the number of genes passing different CV thresholds.

4. Export the Processed Expression Matrix:

The final filtered expression matrix is transposed and saved for WGCNA.

### Usage
Run `01_datExpr_preprocessing.ipynb` to preprocess the expression data. Ensure the input file data/merged_tpm_short_name.csv is available in the specified directory.

## WGCNA Analysis
The WGCNA script 02_WGCNA.R performs the following steps:

1. Install and Load Necessary Packages:

Installs and loads R packages such as WGCNA, here, and matrixStats.

2. Load and Prepare Data:

Reads the preprocessed expression matrix.
Checks for good samples and genes.

3. Pick Soft Threshold:

Determines the appropriate soft threshold power for network construction.

4. Construct Network and Identify Modules:

Constructs the gene co-expression network.
Identifies gene modules.

5. Relate Modules to External Traits:

Associates gene modules with external traits or conditions.

### Usage
Run the R script `02_WGCNA.R` to perform WGCNA. Ensure the preprocessed expression matrix is available in the specified directory.

## References
Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9(1), 559.
