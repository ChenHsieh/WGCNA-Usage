
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("WGCNA", update = TRUE, ask = FALSE)
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  install.packages("matrixStats")
}

library(WGCNA)
library(here)
library(dplyr)
library(matrixStats)

input_file <- here("data", "dynamic_expression_matrix_TPM5_CV30_stem.csv")
soft_threshold_power <- 14
output_directory <- here("output")
max_block_size <- 6000
min_module_size <- 30
merge_cut_height <- 0.25

# Check if the directory exists, if not, create it
if (!file.exists(output_directory)) {
  dir.create(output_directory)
}

datExpr <- read.csv(input_file, row.names = 1)
# datExpr <- as.data.frame(t(datExpr))  # Transpose if necessary to have genes as columns
head(datExpr)

gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # Optionally remove the offending genes and samples from the data
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}


# Function to pick soft threshold
pick_soft_threshold <- function(datExpr) {
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  # Plot the results
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = powers, cex = cex1, col = "red")
  
  return(sft)
}

sft <- pick_soft_threshold(datExpr)


# Function to perform WGCNA
perform_wgcna <- function(datExpr, power, max_block_size, min_module_size, merge_cut_height) {
  net <- blockwiseModules(datExpr, power = power, maxBlockSize = max_block_size,
                          TOMType = "unsigned", minModuleSize = min_module_size,
                          reassignThreshold = 0, mergeCutHeight = merge_cut_height,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE, saveTOMFileBase = "TOM",
                          verbose = 3)
  
  # Convert labels to colors for plotting
  mergedColors <- labels2colors(net$colors)
  num_modules <- length(unique(mergedColors))
  cat("Number of modules detected:", num_modules, "
")
  
  return(list(net = net, mergedColors = mergedColors))
}

wgcna_results <- perform_wgcna(datExpr, soft_threshold_power, max_block_size, min_module_size, merge_cut_height)
net = wgcna_results$net
mergedColors = wgcna_results$mergedColors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# connectivity table output (K)
adjacency = adjacency(datExpr, power = soft_threshold_power)
write.csv(adjacency, file = here(output_directory, "adj_matrix.csv"), row.names = FALSE)

connectivity = intramodularConnectivity(adjacency, net$colors)

connectivity_table = data.frame(
  Gene = colnames(datExpr),
  Module = labels2colors(net$colors),
  Kin = connectivity$kWithin,
  Kout = connectivity$kTotal - connectivity$kWithin,
  Ktotal = connectivity$kTotal
)

head(connectivity_table)

write.csv(connectivity_table, file = here(output_directory, "connectivity_table.csv"), row.names = FALSE)

# net work visualization
MEs = net$MEs
plotEigengeneNetworks(MEs, "Eigengene Network", marDendro = c(0,4,2,0), 
                      marHeatmap = c(3,4,2,2), cex.lab = 0.8, xLabelsAngle = 90)

## module gene number count
# Extract module colors
module_colors <- labels2colors(net$colors)

# Count the number of genes in each module
module_counts <- table(module_colors)

# Create a dataframe
module_counts_df <- as.data.frame(module_counts)
colnames(module_counts_df) <- c("Module", "GeneCount")

# Print the dataframe
print(module_counts_df)

# Optionally, export the dataframe to a CSV file
write.csv(module_counts_df, file = here(output_directory, "module_gene_counts.csv"), row.names = FALSE)

# GOI placement
# Create a dataframe of gene to module mapping
gene_module_mapping <- data.frame(Gene = colnames(datExpr), Module = module_colors)

# get the modules for specific genes
genes_of_interest <- c("PtXaTreH.14G131700", "PtXaTreH.10G125100", "PtXaTreH.06G010600",
                       "PtXaTreH.05G120200", "PtXaTreH.03G064200", "PtXaTreH.02G086600")

module_assignments <- gene_module_mapping[gene_module_mapping$Gene %in% genes_of_interest, ]
print(module_assignments)



# Save results
output_file <- here(output_directory, "wgcna_results.rds")
saveRDS(wgcna_results, output_file)

cat("Results saved to", output_file, "
")



