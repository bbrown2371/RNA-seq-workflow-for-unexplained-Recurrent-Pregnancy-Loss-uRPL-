######################
#### Introduction ####
######################

# This script includes all processes used during my Honors project - 
# Utilising new and improved RNA sequencing (RNA-seq) analysis methods to re-analyse
# existing RNA-seq datasets in unexplained recurrent pregnancy loss (uRPL) research. 

# The code contained here does contain sections where only an example of how 
# an analysis was completed or how a figure can be generated are provided. 
# This has been done to prevent excessive repetition of the same process across
# eleven datasets. 

# There are components of this script that required adjustment
# When a line needed to be adjusted it will have...
# ***CHANGE***
# ...above it

# Install then load core libraries using either: 
# BiocManager (https://bioconductor.org/install) or install.packages
# The tool versions used during our analyses are listed below

# General tools for dataset preparation and manipulation
library(magrittr)          # (version 2.0.3) - For pipe (%>%) and other functional programming tools
library(dplyr)             # (version 1.1.4) - Data manipulation
library(tidyverse)         # (version 2.0.0) - Meta-package including ggplot2, dplyr, tidyr, etc.
library(ggplot2)           # (version 3.5.2) - General plotting 
library(farver)            # (version 2.1.2) - Color space manipulation (used internally by ggplot2 etc.)

# Specialised tools for dataset analysis and visualisation
library(DESeq2)            # (version 1.44.0) - RNA-seq differential expression analysis - used in Sections 2.1 and 3
library(PCAtools)          # (version 2.16.0) - Principal Component Analysis (PCA) and plotting - used in Sections 2.1-2.3
library(clusterProfiler)   # (version 4.12.6) - Functional Annotation of KEGG pathways - used in Section 4
library(pathview)          # (version 1.44.0) - KEGG pathway visualisation - used in Section 4

# The following values are used throughout the script to ensure consistency
# of names and colours assigned to each dataset, tissue and condition
# ***CHANGE***

# Define the names of each dataset and assign a colour for each
datasets <- c("HR", "RF", "RR", "SC", "EBM2", "FIG", "MEFS", "RS", "CTR", "EBM1", "FIP")
dataset_colours <- c("EBM1"="#EA4D96", "FIP"="#BD3BE8", "CTR"="#A99014", "EBM2"="#400D72", "FIG"="#1ABBE0", "HR"="#A8D4F5",
                     "MEFS"="#18CD61", "RF"="#841063", "RR"="#55BB16", "RS"="#E4EE72", "SC"="#EC8A5F")
dataset_shapes <- c('EBM1'=0, 'FIP'=13, 'CTR'=9, 'EBM2'=4, 'FIG'=2, 'HR'=3, 'MEFS'=15, 'RF'=6, 'RR'=7, 'RS'=8, 'SC'=18)

# Define colours associated with each reproductive tissue type
tissue_colours <- c("endometrial"="#96F3D2", "decidual"="#C8BBF7", "chorionic villus"="#92F184", "endometrial gland"="#F9CDF9")

# Define the shapes and colours to represent the experimental condition/group
condition_shapes <- c("control"=16, "uRPL"=17)
condition_colours <- c("control"="blue", "uRPL"="brown2")

################################################
#### Section 1: Initial Dataset Preparation ####
################################################

# Load in Salmon gene-level counts files generated using the nf-core/rnaseq pipeline (https://nf-co.re/rnaseq/3.18.0/) 
# and associated sample metadata 
# Ensure that the counts matrix and metadata:
# - Include matching sample identifiers
# - Are in the same order
# - Have metadata variables as either factors or numeric types (required for DESeq2)

# In the code used below the following inputs as named:
# Gene-level counts matrix (containing all samples) = "counts"
# Sample metadata (containing all datasets) = "samples_table"

#######################################################
#### Section 2: Principal Component Analysis (PCA) ####
#######################################################

### --- Combined dataset PCA ---

# Create DESeq2 dataset from full counts matrix and sample metadata and run DESeq2 pipeline
# Note: condition = experimental group
dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts, colData = samples_table, design = ~condition)
deseq <- DESeq(dds)

# Normalise the count data by applying a Variance Stabilising Transformation (vst)
vst <- assay(vst(deseq))

# Perform PCA on the vst-transformed data 
p <- pca(vst, metadata = samples_table)

# Generate a data frame that contains the rotated PCA coordinates (principal component scores)
pc_df <- p$rotated
# Add a variable to the data frame that labels specified data points (i.e. outliers)
pc_df$label <- ifelse(rownames(pc_df) %in% c("EBM1_C3", "MEFS_C1", "RS_C3"), rownames(pc_df), NA)
# Extract % variance explained by PC1 and PC2 for labeling axes (rounded to one decimal place but can be changed)
pc1_var <- round(p$variance["PC1"], 1)
pc2_var <- round(p$variance["PC2"], 1)

# Create PCA plot
ggplot(pc_df, aes(x = PC1, y = PC2),    ) + theme_light(base_size = 20) + 
  # Point geometry coloured by dataset and shaped by experimental condition
  geom_point(aes(color = p$metadata$Dataset, shape = p$metadata$condition), size = 3) +
  # Colour scale for dataset
  scale_color_manual(values = dataset_colours) + 
  # Add a legend for conditions
  scale_shape_manual(values = condition_shapes, name = "Condition") +
  # Add labels to selected points using variable designed above
  geom_text_repel(aes(label = label), size = 4, na.rm = TRUE, point.padding = 0.3, box.padding = 0.7) +  
  # Axis labels with % variance 
  labs(color = "Dataset", shape = "Group", 
       x = paste0("PC1 (", pc1_var, "% variance)"), 
       y = paste0("PC2 (", pc2_var, "% variance)")) +
  # Allow the use of second colour scale (used to add tissue type ellipses)
  ggnewscale::new_scale_color() + 
  # Add 95% confidence intervals for each tissue type
  stat_ellipse(aes(group = p$metadata$tissue_type, fill = p$metadata$tissue_type, colour = p$metadata$tissue_type), 
               size = 1, geom = "polygon", alpha = 0.15, level = 0.95) + 
  # Define the fill and border colour for tissue type ellipses
  scale_fill_manual(values = tissue_colours, name = "Tissue Type") + 
  stat_ellipse(aes(group = p$metadata$tissue_type, color = p$metadata$tissue_type), 
               geom = "path", size = 0.5, level = 0.95) +
  scale_fill_manual(values = tissue_colours, name = "Tissue Type") +
  # Hides duplicate legend for the confidence interval border colours
  scale_color_manual(values = tissue_colours, guide = "none") 

### --- PCA by tissue type subset (repeat PCA on tissue-specific subsets) ---
## -- endometrial samples --
endo_samples <- samples_table %>% filter(Dataset %in% c('SC','RF','RR','HR'))
endo_counts <- counts %>% select(matches("SC|RF|RR|HR"))

dds_endo <- DESeq2::DESeqDataSetFromMatrix(countData=endo_counts, colData = endo_samples, design = ~ condition)
deseq_endo <- DESeq(dds_endo)
vst_endo <- assay(vst(deseq_endo))
p_endo <- pca(vst_endo, metadata = endo_samples)

## -- Decidual samples --
deci_samples <- samples_table %>% filter(Dataset %in% c('RS','FIG','MEFS','EBM2'))
deci_counts <- counts %>% select(matches("RS|FIG|MEFS|EBM2"))

dds_deci <- DESeq2::DESeqDataSetFromMatrix(countData=deci_counts, colData = deci_samples, design = ~ condition)
deseq_deci <- DESeq(dds_deci)
vst_deci <- assay(vst(deseq_deci))
p_deci <- pca(vst_deci, metadata = deci_samples)

## -- Chorionic villus samples --
villus_samples <- samples_table %>% filter(Dataset %in% c('EBM1','FIP','CTR'))
villus_counts <- counts %>% select(matches("EBM1|FIP|CTR"))

dds_villus <- DESeq2::DESeqDataSetFromMatrix(countData=villus_counts, colData = villus_samples, design = ~ condition)
deseq_villus <- DESeq(dds_villus)
vst_villus <- assay(vst(deseq_villus))
p_villus <- pca(vst_villus, metadata = villus_samples)

# At this point, individual PCA plots were generated for each tissue type subset
# using the ggplot structure used above with <p> substituted for <p_endo>, <p_deci>, and <p_villus>.

#############################################
#### Section 2.1: Generating PCA objects ####
#############################################

### --- PCA of each dataset ---

# Initialise empty lists to store the DESeq2 and PCA objects per each dataset
ds <- list()
p_list <- list()

# Generate DESeq2 and PCA objects for each dataset
for (name in datasets) {
  # Subsets metadata and counts matrix to only include samples from the current dataset
  samples <- samples_table %>% filter(Dataset %in% c(name))
  counts_sub <- counts %>% select(matches(name))
  
  # Creates a DESeq object using the raw counts and metadata and runs DESeq2 pipeline
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts_sub, colData = samples, design = ~ condition)
  deseq <- DESeq(dds, quiet = TRUE)
  ds[[name]] <- deseq  # Save the DESeq2 object 
  
  # Normalises the count data by applying a Variance Stabilising Transformation (vst)
  vst <- assay(vst(deseq))
  
  # Performs PCA on the vst-transformed data 
  p <- pca(vst, metadata = samples)
  p_list[[name]] <- p # Save the PCA object
}

####################################################
#### Section 2.2: Principal Component Retention ####
####################################################

# The results from each method can be used to determine a consensus principal component (PC) 
# to retain per dataset, which can then be used to identify confounding factors

# Initialize a list to store retained PCs for each method
PCs <- list()

# Run the loop to calculate the retained PCs for each dataset using six methods
for (name in names(ds)) {
  # Retrieves the DESeq2 object from the named list (ds) and isolates the associated metadata
  data <- ds[[name]] 
  samples <- samples_table %>% filter(Dataset %in% c(name))
  
  vst <- assay(vst(data))
  p <- pca(vst, metadata = samples)
  
  # 1. (Horn's) Parallel Analysis
  # Retains PCs that explain more variance than random noise
  horn_PC <- parallelPCA(vst, threshold = 0.05)
  PCs$horn[[name]] <- horn_PC$n
  
  # 2. Scree test
  # Finds elbow point in variance curve
  elbow_PC <- findElbowPoint(p$variance)
  PCs$elbow[[name]] <- elbow_PC
  
  # Calculates eigenvalues = variance explained * number of samples
  eigenvalues <- (p$variance/100) * nrow(p$loadings) 
  
  # 3. Jolliffe's modification of Kaiser-Guttman 
  # Retains PCs with eigenvalues greater than 0.7 * mean eigenvalue
  jkg_PC <- which(eigenvalues > (mean(eigenvalues)* 0.7))
  PCs$jkg[[name]] <- jkg_PC
  
  # 4 . Kaiser-Guttman method
  # Retains PCs with eigenvalues greater than the mean eigenvalue
  kg_PC <- which(eigenvalues > mean(eigenvalues))
  PCs$kg[[name]] <- kg_PC
  
  # 5. Cumulative variance
  # Retains PCs explaining more than 80% of the total variance
  PCs$cumsum80[[name]] <- which(cumsum(p$variance) > 80) [1]
  
  # 6. Modified Broken stick model
  # Compares PC group variances against broken stick thresholds to determine which PCs to retain
  N <- nrow(p$metadata) # Number of samples in the dataset
  threshold <- 1.5 # Note: value is somewhat arbitrary but recommended to use a value between 1 and 2
  
  # Estimates the sampling error (δλ) for each eigenvalue
  delta_lambda <- eigenvalues * sqrt(2 / N) 
  
  # Groups eigenvalues into subspaces based on proximity (i.e. when spacing ≤ threshold) 
  groups <- list()
  current_group <- c(1)
  for (i in 2:length(eigenvalues)) {
    spacing <- abs(eigenvalues[i - 1] - eigenvalues[i])
    if (spacing <= threshold) {
      current_group <- c(current_group, i)
    } else {
      groups[[length(groups) + 1]] <- current_group
      current_group <- c(i) }   }
  groups[[length(groups) + 1]] <- current_group
  
  # Calculates the total variance explained by each group of eigenvalues (subspace)
  Wi_sums <- sapply(groups, function(g) sum(eigenvalues[g]))
  
  # Computes Broken Stick expectations for each subspace
  mod_bs <- function(n) {
    sapply(1:n, function(k) sum(1 / (k:n)) / n)   }
  bs_values <- mod_bs(length(Wi_sums)) * sum(eigenvalues)   
  
  PCs$mod_bs[[name]] <- bs_values   }

# Example usage: View retained PCs using PCs$<dataset_name>

########################################################
#### Section 2.3: Confounding Factor Identification ####
########################################################

# Note: this section shows an example of the methods and code used to identify confounding variables
# This code is generalised and was adjusted to investigate each dataset and/or PC as needed

## Option 1: Eigencor plots - to be used when they are at least two numeric variables and 
## more than four samples per group
# ***CHANGE***
eigencorplot(p_list$<dataset_name>, # Replace with actual dataset name e.g. p_list$HR
             metavars = c("age", "BMI", "prev_losses", "live_births", "gestational_age", "cycle_length"), 
             col = c("darkblue", "blue", "black", "red", "darkred"), colCorval = "white", scale = TRUE, 
             main = "Eigencor plot")
# metavars = list of numeric variables to check for correlation with PCs (adjust as needed)

## Option 2: Manual Correlation - to be used if a dataset has only one numeric variable or
## has a small sample size (less than 4 samples per group)
# Extract PC scores for the dataset
# ***CHANGE***
pc_df <- p_list$<dataset_name>$rotated
# Example usage: Check the correlation between PC1 and age
cor.test(pc_df$PC1, p_list$<dataset_name>$metadata$age)

## Option 3: Colouring a PCA plot according to a specified variable
# This can be used when you want to investigate the relationship between a factor variable 
# and gene expression (e.g. determine if cycle day at biopsy is affecting the gene expression profiles)
# Note: the example below shows how to generate a PCA plot with each sample coloured according to cycle day at biopsy
# for dataset SC, which was used to generate Supplementary Figure S3

# Extract PC scores  and variance explained by PC1 and PC2 for the dataset
# ***CHANGE***
pc_df <- p_list$<dataset_name>$rotated
meta <- p_list[[<dataset_name>]]$metadata
pc1_var <- round(p_list$<dataset_name>$variance["PC1"], 1)
pc2_var <- round(p_list$<dataset_name>$variance["PC2"], 1)

# Generate a PCA plot coloured by cycle day at biopsy (not used for my analyses)
ggplot(pc_df, aes(x = PC1, y = PC2)) + theme_light(base_size = 20) + 
  geom_point(aes(colour = meta$cycle_day, shape = meta$condition), size = 2.5) + 
  scale_shape_manual(values = c("control"=16, "uRPL"=17), name = "Condition") +
  stat_ellipse(aes(group = meta$cycle_day, fill = meta$cycle_day, colour = meta$cycle_day), 
               size = 0.7, geom = "polygon", alpha = 0.1, level =0.95) +
  scale_color_manual(values = c("LH+8"="cornflowerblue", "LH+10"="brown2", "LH+7"="limegreen", "LH+6"="goldenrod", "LH+9"="orchid")) + 
  scale_fill_manual(values = c("LH+8"="cornflowerblue", "LH+10"="brown2", "LH+7"="limegreen", "LH+6"="goldenrod", "LH+9"="orchid"), name = "Reported Cycle Day") +
  labs(color = "Reported Cycle Day", shape="Condition", x = paste0("PC1 (", pc1_var, "% variance)"), y = paste0("PC2 (", pc2_var, "% variance)"))


###########################################################################
#### Section 2.4: PCA plots with arrows representing numeric variables ####
###########################################################################

# The arrows generated represent the *direction* of correlation between numeric metadata 
# variables (e.g., age, BMI) and the principal components. 
# Note: arrows only show *direction*, NOT statistical significance or strength of the correlation

# Define which numeric variables are to be plotted as directionality arrows for each dataset 
# Note: variables must be numeric and have values for all samples otherwise this will throw an error
# ***CHANGE (AS NEEDED)***
# e.g. the confounding factors identified above can be added to the DESeq2 design
metadata_vars_list <- list(HR = c("prev_losses", "age", "live_births"), RF = NULL,  RR = NULL,  SC = c("prev_losses"),
                           EBM2 = NULL,  FIG = NULL,  MEFS = NULL,  RS = c("live_births"),
                           CTR = NULL,  EBM1 = NULL,  FIP = NULL)

# Generate the function that generates a PCA plot with directionality arrows for 
# the specified numeric metadata variables per dataset
plot_pca_with_arrows <- function(datasets, metadata_vars, p_list) {
  # Generates a data frame that contains the PCA coordinates and one containing the associated metadata
  pc_df <- p_list[[datasets]]$rotated
  meta <- p_list[[datasets]]$metadata
  # Extracts % variance explained by PC1 and PC2 for labeling axes
  pc1_var <- round(p_list[[datasets]]$variance["PC1"], 1)
  pc2_var <- round(p_list[[datasets]]$variance["PC2"], 1)
  
  # Creates the base PCA plot
  p_base <- ggplot(pc_df, aes(x = PC1, y = PC2)) + theme_light(base_size = 20) + 
    # Point geometry coloured by dataset and shaped by experimental condition
    geom_point(aes(color = meta$condition, shape = meta$condition), size = 4) + 
    # Add 95% confidence intervals for each condition
    stat_ellipse(aes(group = meta$condition, fill = meta$condition, colour = meta$condition), 
                 size = 0.7, geom = "polygon", alpha = 0.1, level = 0.95) + 
    # Define the colour and shape scales 
    scale_shape_manual(values = condition_shapes, name = "Condition") +
    scale_color_manual(values = condition_colours, name = "Condition") + 
    scale_fill_manual(values = condition_colours, name = "Condition") +
    # Axis labels with % variance 
    labs(x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"), title = datasets)
  
  # Checks if the dataset has numeric variables to plot as arrows
  # Arrows will only be added is metadata_vars is not empty for the current dataset
  if (!is.null(metadata_vars) && length(metadata_vars) > 0) {
    arrow_data <- data.frame(Variable = metadata_vars, PC1 = NA, PC2 = NA)
    
    # Loops through each variable and compute the Pearson correlation coefficient with PC1 and PC2
    for (var in metadata_vars) {
      arrow_data[arrow_data$Variable == var, "PC1"] <- cor(meta[[var]], pc_df$PC1, use = "complete.obs")
      arrow_data[arrow_data$Variable == var, "PC2"] <- cor(meta[[var]], pc_df$PC2, use = "complete.obs") }
    
    # Scales the arrow vectors (to make them more visible on the PCA plot)
    # Note: this scales arrow length uniformly — it does not preserve magnitude.
    # ***CHANGE (AS NEEDED)***
    arrow_data$PC1 <- arrow_data$PC1 * 50
    arrow_data$PC2 <- arrow_data$PC2 * 50
    
    # Adds variable loading arrows to base PCA plot
    p_base <- p_base + 
      # Draws arrows from origin (0,0) to the scaled PC1/PC2 coordinates
      geom_segment(data = arrow_data, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                   arrow = arrow(length = unit(0.3, "cm")), color = "black") +
      # Adds labels to the tip of each arrow
      geom_text(data = arrow_data, aes(x = PC1, y = PC2, label = Variable), vjust = -0.5, hjust = 0.75, size = 6)   }
  
  return(p_base)   }

# Loop the function through each dataset and save the PCA plots in a list 
pca_arrow_plots <- lapply(datasets, function(ds) {
  plot_pca_with_arrows(ds, metadata_vars_list[[ds]], p_list)  })

# Rename each plot using the named list for clarity
names(pca_arrow_plots) <- datasets

# Example usage: View PCA plot for a dataset using pca_arrow_plots[["<dataset_name>"]]

#############################################################
#### Section 3: Differential Expression Analysis (DEA)   ####
#############################################################

# Note: dataset HR was removed from all analyses starting here

## Complete DEA on each dataset using a function containing the DESeq2 pipeline

# Create a named list containing the corresponding DESeq2 design for each dataset 
# including previously identified confounding factors
# ***CHANGE (AS NEEDED)***
designs <- list( RF = ~ condition, RR = ~ condition, SC = ~ condition,
                 EBM2 = ~ condition, FIG = ~ condition, MEFS = ~ condition, RS = ~ live_births + condition,
                 CTR = ~ condition, EBM1 = ~ condition, FIP = ~ condition )

# Generate a function that will run each dataset through the DESeq2 pipeline
run_deseq <- function(dataset, design_formula, counts, samples_table) {
  # Subsets metadata and count data for the current dataset
  dataset_samples <- samples_table %>% filter(Dataset %in% c(name))
  dataset_counts <- counts %>% select(matches(dataset))
  
  # Creates DESeq2 object and runs the DESeq2 pipeline
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = dataset_counts, colData = dataset_samples, design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  
  return(dds)  }

# Run the DESeq2 function on each dataset listed in the design element aobve
# and store results in a named list (DEA)
DEA <- lapply(names(designs), function(ds) {
  run_deseq(dataset = ds, design_formula = designs[[ds]], 
            counts = counts, samples_table = samples_table)  })

# Rename the results to match each dataset name to improve clarity
names(DEA) <- names(designs)

## Summarise the differentially expressed genes (DEGs) for each dataset

# Create an empty list to store DEG results per dataset
deg_up_lists <- list()
deg_down_lists <- list()
deg_lists <- list()

# Identify the significant (FDR-adjusted p-value < 0.05) DEGs for each dataset
for (name in names(DEA)) {
  dds_data <- DEA[[name]] # Gets the DESeq2 object
  results <- results(dds_data, alpha = 0.05)
  # Filters the results for significant DEGs (FDR-adjusted p-value < 0.05)
  total_degs <- rownames(subset(results,  padj < 0.05))
  deg_lists[[name]] <- total_degs
  
  # Filters the results for significant DEGs and splits by direction of differential expression
  # (FDR-adjusted p-value < 0.05 and positive or negative log2 fold change) and extract the gene names only
  sigDEG_up <- subset(results,  padj < 0.05 & log2FoldChange > 0)
  sig_genes_up <- rownames(sigDEG_up)
  deg_up_lists[[name]] <- sig_genes_up
  
  sigDEG_down <- subset(results,  padj < 0.05 & log2FoldChange < 0)
  sig_genes_down <- rownames(sigDEG_down)
  deg_down_lists[[name]] <- sig_genes_down   }

# Combine the DEGs from each study into a data frame
summary_DEG_df <- do.call(rbind, lapply(names(deg_up_lists), function(name) {
  # Gets the DEG list and collapses it into one string (separated by commas) for each dataset
  up_DEGs <- deg_up_lists[[name]]
  up_DEG_string <- paste(up_DEGs, collapse = ", ")
  up_total <- length(up_DEGs)
  down_DEGs <- deg_down_lists[[name]]
  down_DEG_string <- paste(down_DEGs, collapse = ", ")
  down_total <- length(down_DEGs)
  data.frame(Dataset = name, Total_Upregulated = up_total, Upregulated_DEGs = up_DEG_string, 
             Total_Downregulated = down_total, Downregulated_DEGs = down_DEG_string, stringsAsFactors = FALSE)    }))

# Combine both up and down DEG strings into one vector
all_DEG_strings <- c(summary_DEG_df$Upregulated_DEGs, summary_DEG_df$Downregulated_DEGs)
# Split the strings into individual gene names
all_DEGs <- unlist(strsplit(all_DEG_strings, split = ",\\s*"))
# Count the frequency of each DEG across all datasets
DEG_freq_table <- as.data.frame(sort(table(all_DEGs), decreasing = TRUE))
# Rename columns for clarity
colnames(DEG_freq_table) <- c("DEG", "Frequency")

## Visualise the overlap of significant DEGs across each tissue type using Venn diagrams

# Load library for Venn diagram plotting 
library(ggvenn)

# Group datasets by tissue type into named lists containing the list of DEGs for each dataset
endo_DEGs <- list(SC = deg_lists$SC, RR = deg_lists$RR, RF = deg_lists$RF)
decidua_DEGs <- list(EBM2 = deg_lists$EBM2,RS = deg_lists$RS, FIG = deg_lists$FIG, MEFS = deg_lists$MEFS)
villus_DEGs <- list(EBM1 = deg_lists$EBM1, CTR = deg_lists$CTR, FIP = deg_lists$FIP)

# Create Venn diagrams for each tissue type
# Other tissue types can be visualised by replacing 'endo_DEGs' with 'decidua_DEGs' or 'villus_DEGs'
ggvenn(endo_DEGs, fill_color = dataset_colours) 

## Visualise the normalised read counts of a gene using box plots

# Specify the gene to be plotted
# ***CHANGE (AS NEEDED)***
gene <- "CDKN2B"

# Create a list to save the ggplot box plots for each dataset
boxplot_list <- list()

# Loop over each dataset in the named list 'DEA' to generate box plots facetted by dataset
for (name in names(DEA)) {
  # Extracts the DESeq2 object for the current dataset
  ds <- DEA[[name]]
  
  # Checks if the gene is found in the row names, 
  # if not (ie. not in the dataset) the dataset will be skipped
  if (!gene %in% rownames(ds)) {
    message(paste("Skipping", name, "- gene not found"))
    next  }
  
  # Extracts the normalised read counts for the dataset
  norm_counts <- counts(ds, normalized = TRUE)
  
  # Gets the DESeq2 results and extract the adjusted p-value (padj) and 
  # log2 fold change (lfc) for the gene
  res <- results(ds)
  padj <- signif(res[gene, "padj"], 3)
  lfc <- signif(res[gene, "log2FoldChange"], 3)
  
  # Prepares a data frame containing the normalised read counts 
  # for the gene and the sample metadata
  plot_df <- norm_counts[gene, , drop=FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
    left_join(as.data.frame(colData(ds)) %>% rownames_to_column("sample"), by = "sample")
  
  # Creates the ggplot box plot showing outliers
  p <- ggplot(plot_df, aes(x = condition, y = expression, fill = condition)) + theme_light() + 
    geom_boxplot(outliers = FALSE) + geom_jitter(width = 0.2, size = 1, alpha = 0.6) + 
    scale_fill_manual(values = c("control" = "cornflowerblue", "uRPL" = "brown2")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none") +
    labs(title = paste("Dataset", name), y = "Normalised Counts", subtitle = paste("log2FC:", lfc, "| padj:", padj))
  
  boxplot_list[[name]] <- p   }

# Example usage: box plots can be visualised individually using boxplot_list[[<dataset_name>]] 
# alternatively all plots can be viewed together using patchwork::wrap_plots(boxplot_list)

############################################################################
#### Section 4: Functional Annotation: KEGG pathway enrichment analysis ####
############################################################################

# Load "org.Hs.eg.db" package to enable mapping of Entrez Gene IDs 
library(org.Hs.eg.db)

# Create empty lists to store the KEGG results and diagrams for each study
kegg_plot <- list()
kegg_table <- list()
DEA_res <- list()

## Complete enrichment analysis on each dataset
# Loop through for each dataset in the named list 'DEA' to generate
# KEGG results tables and dot plots
for (name in names(DEA)) {
  # Extracts the DESeq2 object and filters for significant DEGs 
  ds <- DEA[[name]]
  res <- results(ds, alpha = 0.05)
  resgenes <- subset(res, padj < 0.05)
  DEA_res[[name]] <- resgenes
  
  # Converts gene symbols to Entrez IDs using the org.Hs.eg.db database
  # Any datsets with no significant DEGs are skipped
  if (nrow(resgenes) == 0) next
  gene_list <- bitr(rownames(resgenes), fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  
  # Performs KEGG enrichment analysis and generates a data frame of the results 
  # Any datasets where no valid gene mappings were found are skipped
  if (nrow(gene_list) == 0) next
  kk <- enrichKEGG(gene = gene_list$ENTREZID,organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  kktable <- data.frame(kk)
  
  # Generates a dot plot of the top enriched pathways
  # Any datasets where KEGG enrichment analysis was unsuccessful 
  # (i.e. has no significantly enriched pathways) are skipped
  if (!is.null(kk) && nrow(kk) > 0) {
    kegg_dot <- dotplot(kk, showCategory = 14, font.size = 12) + ggtitle(paste("Dataset", name))
    kegg_plot[[name]] <- kegg_dot
    kegg_table[[name]] <- kktable  }  }

# Remove empty plots from the plot list (in case some datasets were skipped or failed)
kegg_plot <- kegg_plot[!sapply(kegg_plot, is.null)]

# Example usage: plots can be visualised together using patchwork::wrap_plots(kegg_plot)

## Summarise the results of the KEGG enrichment analysis

# Combine the KEGG results from each study into a single data frame
summary_KEGG_df <- do.call(rbind, lapply(names(kegg_table), function(name) {
  # Gets the KEGG pathway list and collapses it into one string (separated by commas)
  # for each dataset
  Pathways <- kegg_table[[name]]$Description
  Pathway_string <- paste(Pathways, collapse = ", ")
  Total <- length(Pathways)
  data.frame(Dataset = name, Count = Total, Pathways = Pathway_string, stringsAsFactors = FALSE)    }))

# Split the pathway strings and count the frequency of each pathway across all studies
all_pathways <- unlist(strsplit(summary_KEGG_df$Pathways, split = ",\\s*"))
# Create and order the table by decreasing frequency
pathway_freq_table <- as.data.frame(sort(table(all_pathways), decreasing = TRUE))
# Rename the columns for clarity
colnames(pathway_freq_table) <- c("Pathway", "Frequency") 

#### Generate a mirrored bar chart to summarise DEA within a chosen KEGG pathway ####

# Upload a file containing the genes within the chosen KEGG pathway 
# ***CHANGE***
pathway_genes <- read.csv("./Cytokinepathway_genes.csv", header=TRUE, stringsAsFactors = FALSE)

dea_pathway <- list()

# Loop through each dataset that generated KEGG enrichment results above
for (name in names(kegg_table)) {
  
  # Creates a data frame from the DESeq results for the specified dataset
  dea_subset <- as.data.frame(results(DEA[[name]]))
  # Isolates the genes existing with the chosen KEGG pathway within 
  pathway_subset <- dea_subset %>% filter(rownames(dea_subset) %in% pathway_genes$gene)
  
  dea_pathway[[name]] <- pathway_subset     }

# Combine the result subsets into one data frame adding an id column according to dataset
all_res <- bind_rows(dea_pathway, .id = "dataset")

# Define the DEG significance cut-off (i.e. the adjusted p-value)
# ***CHANGE (AS NEEDED)***
lfc_cutoff <- 2

# Generate a data frame using the combined results that adds a column
# specifying up or down regulation based on the log2FoldChange
plot_data <- all_res %>% filter(!is.na(padj)) %>%
  # Labels each gene as either above or below the log2FoldChange cut-off
  mutate(lfc = ifelse(abs(log2FoldChange) > lfc_cutoff, "|log2FC| > 2", "|log2FC| < 2"),
         # Classifies genes as up (positive log2FoldChange) or down (negative log2FoldChange) regulated
         Regulation = ifelse(log2FoldChange > 0, "Up", "Down"),
         # Bins adjusted p-values (padj) into significance levels and labels each level
         # i.e. -Inf to 0.001 label is '<= 0.001' and the 0.001 to 0.01 label is '0.001-0.01'
         sig_bin = cut(padj,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, 1.00),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", "n.s")),
         # Defines the variable ('x') which indicates the direction of DE for each gene
         # +1 for up-regulation and -1 for down-regulation
         x = ifelse(Regulation == "Down", -1, 1))

# Determine the number of up and down-regulated DEGs for each dataset (optional)
deg_counts <- plot_data %>% group_by(dataset, Regulation, lfc, sig_bin) %>%
  summarise(Count = n(), .groups = "drop")

# Generate a mirrored bar plot using ggplot2 with each dataset separated on y-axis
ggplot(plot_data, aes(x = x, y = dataset, fill = sig_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  # Adds a vertical line at the origin to separate the up and down-regulated genes (optional)
  geom_vline(xintercept = 0, linewidth = 0.6) +
  # Specifies the colours for each significance level 
  scale_fill_manual(name = "Adj. p-value", values = c("<0.001" = "#453781FF", "0.001-0.01" = "#238A8DFF", 
                                                      "0.01-0.05" = "#55C667FF", "n.s" = "#FDE725FF")) +
  labs(x = "Number of DEGs", y = "Dataset") + theme_light(base_size = 13)


#### Visualise the cytokine-cytokine receptor pathway on a KEGG map ####

# Load "org.Hs.eg.db" package to enable mapping of Entrez Gene IDs 
library(org.Hs.eg.db)

# Specify which KEGG pathway to be mapped using the KEGG pathway entry ID
# ***CHANGE*** 
kegg_pathway <- "hsa04060"

pathview_datasets <- list()
KEGG_maps <- list()

# Loop through each dataset specified in the kegg_table list to generate a dataframe of log2FoldChange
# and a KEGG pathway visualisation showing the regulation of each gene within the specified pathway
# Note: kegg_table will only include the datasets where KEGG enrichment results were generated
for (name in names(kegg_table)){
  # Extract log2FoldChange results for specified dataset from DEA results
  dataset <- as.data.frame(results(DEA[[name]]))
  dataset_lfc <- dataset[, c("log2FoldChange"), drop = FALSE]
  
  # Create a column containing the gene symbols (i.e. the row names of the DEA results)
  dataset_lfc$SYMBOL <- rownames(dataset)
  
  # Map the gene symbols to Entrez IDs 
  dataset_lfc$ENTREZID <- mapIds(org.Hs.eg.db, keys = dataset_lfc$SYMBOL, column = "ENTREZID", 
                                 keytype = "SYMBOL", multiVals = "first")
  
  # Filter out genes with no Entrez ID mapping to rename the rows as Entrez IDs 
  # rather than gene symbols (NAs cannot be used to name rows) 
  dataset_lfc <- dataset_lfc[!is.na(dataset_lfc$ENTREZID), ]
  rownames(dataset_lfc) <- dataset_lfc$ENTREZID
  
  # Remove the Entrez and gene symbol columns so only columns containing LFCs remain
  # with the Entrez IDs present as the row names
  dataset_lfc <- subset(dataset_lfc, select=-c(SYMBOL, ENTREZID))
  
  pathview_datasets[[name]] <- dataset_lfc
  
  # Generate KEGG pathway map for the specified dataset
  # Note: the images are automatically saved to your directory 
  KEGG_map <- pathview(gene.data = dataset_lfc, pathway.id = kegg_pathway, species = "hsa",
                       kegg.native = TRUE, multi.state  = TRUE, same.layer = TRUE, out.suffix=paste(name), 
                       low = list(gene = "orange"), mid = list(gene = "gray"), high = list(gene = "cornflowerblue"))
  
  KEGG_maps[[name]] <- KEGG_map    }

# Note: KEGG pathway maps are automatically saved to working directory 
# Example usage: The information used to generate each plot can be accessed
# using KEGG_maps[[<dataset_name>]] 

#############################################################################################
#### Investigating the impact of dataset variability on differential expression analysis ####
#############################################################################################

## Complete differential expression analysis on the combined dataset (minus dataset HR)

# Define datasets to include (excluding HR)
include_datasets <- c("SC", "RR", "EBM2", "RF", "FIG", "MEFS", "RS", "CTR", "EBM1", "FIP")

# Subset samples and counts
comb_samples <- samples_table %>% filter(Dataset %in% include_datasets)
comb_counts  <- counts %>% select(matches(paste(include_datasets, collapse = "|")))

# Run DESeq2
dds <- DESeq2::DESeqDataSetFromMatrix(countData=comb_counts, colData = comb_samples, design = ~condition)
deseq <- DESeq(dds)
res_all <- results(deseq, alpha = 0.05)

# Helper to extract DEGs
get_DEGs <- function(res, direction = c("up", "down")) {
  direction <- match.arg(direction)
  if (direction == "up") {
    res <- res %>% subset(padj < 0.05 & log2FoldChange > 0)} 
  else {  res <- res %>% subset(padj < 0.05 & log2FoldChange < 0) }
  genes <- rownames(res)
  list(genes = genes, string = paste(genes, collapse = ", "), total = length(genes))   }

# Extract up- and downregulated DEGs
up   <- get_DEGs(res_all, "up")
down <- get_DEGs(res_all, "down")

# Summarise the results 
summary_DEG_df <- data.frame(Total_Upregulated   = up$total, Upregulated_DEGs    = up$string,
                             Total_Downregulated = down$total, Downregulated_DEGs  = down$string, stringsAsFactors    = FALSE)

# Prepare a Venn diagram
library(ggvenn)

# Extract the concordant DEGs identified by comparing the individual datasets
concordant_DEGs <- DEG_freq_table %>% filter(Frequency >= 2) %>% pull(DEG) %>% as.character()


# Extract all significant DEGs from the combined dataset regardless of direction of differential expression
sig_all         <- subset(res_all, padj < 0.05)
sig_all_genes   <- rownames(sig_all)

# Generate a Venn diagram comparing the signficant DEGs identified by each method
comb_DEGs <- list(concordant_DEGs = concordant_DEGs, overall_DEGs = sig_all_genes)
ggvenn(comb_DEGs)

# Identify which DEGs were common to the combined dataset and the concordant DEGs identified 
# by comparing the results from each dataset
intersection <- intersect(concordant_DEGs, sig_all_genes)

## Complete functional annotation of KEGG pathways to identify enriched pathways 

# Load "org.Hs.eg.db" package to enable mapping of Entrez Gene IDs 
library(org.Hs.eg.db)

# Run KEGG pathway enrichment analysis and generate a dot plot
comb_gene_list <- bitr(rownames(comb_sig), fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
comb_kk <- enrichKEGG(gene = comb_gene_list$ENTREZID,organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(comb_kk, showCategory = 14, font.size = 14)

