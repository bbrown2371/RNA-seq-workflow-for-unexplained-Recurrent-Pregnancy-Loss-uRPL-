# Expected Outputs of RNA-seq workflow

## Section 1: Intial Dataset Preparation 
Once the data has been read into R studio, you should have a file containing the Salmon gene expression matrices for each sample (generated from the nf-core/rnaseq pipeline) and a table containing the associated metadata for each sample. 
> [!NOTE]
> Ensure that the counts matrix and metadata:
>  - Include matching sample identifiers (e.g. SC_C1)
>  - Samples are in the same order across both files 
>  - Have all metadata variables set as either factors or numeric variables

## Section 2: Principal Component Analysis
### Section 2.1: Generate PCA objects
This section must be completed before attempting any other Section 2 components.\
[DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and [PCAtools](https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html) are used to generate the objects that are the required input for steps other components of Section 2.

### Section 2.2: Principal Component (PCs) Retention
Generates a named list that contains the retained PCs suggested by six different PC retention methods.\
These can then be used to determine which PC should be used as the limiting PC (i.e. a consensus PC) to identify confounding factors in Section 2.3.\
Depending on the method it may produce a single PC which is the limiting PC or it may return a set of PCs that should be retained according to that method (see below).
| PC retention method | Example Output |
| :------- | :------ |
| Horn's Parallel Analysis | Limiting PC (e.g. 3) |
| Scree test (Elbow point) | Limiting PC (e.g. 3) |
| Jolliffe's modification of the Kaiser-Guttman method | Retained PCs (e.g. 1 2 3 4 5) |
| Kaiser-Guttman method | Retained PCs (e.g. 1 2) |
| Cumulative Variance (>80%) | Limiting PC (e.g 8) |
| Modified Broken Stick Method | Retained PCs (e.g. 1 2 16 17 18 19) |

### Section 2.3: Confounding factor identification 
Using the consensus/limiting PC determined in Section 2.2, [PCAtools](https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html) can generate eigencor plots to determine which numeric variables are confounding factors.\
For example:\
According to this eigencor plot (assuming the limiting PC is the third PC), the only confounding factor is the number of previous losses ("prev_loss").

<img width="630" height="440" alt="SC_eigen" src="https://github.com/user-attachments/assets/f8d96779-b8bd-4610-9d30-1c7b397e29b7" />

When an eigencor plot is not possible, manual Pearson's correlation can be as such:
```
# Note: p_list$,dataset_name> represents the PCA object for a specified dataset
pc_df <- p_list$<dataset_name>$rotated
cor.test(pc_df$PC1, p$metadata$age)
```

To check if non-numeric or factor variables are confounding factor for a dataset, a PCA plot can be generated where each sample is coloured according to that variable.\
For example:
This plot shows each sample of a dataset coloured according to the reported cycle day at tissue sample biopsy.

<img width="600" height="400" alt="SC_ggplot_cycle" src="https://github.com/user-attachments/assets/7b49ffdd-9dbe-42ca-aa34-dd91f8289d91" />

### Section 2.4: Generating PCA plot with arrows for numeric variables 
This process has been designed to allow you to select which numeric clinical variable/s you wish to show on the PCA plot. 

For example:\
These PCA plots show several different clinical variables plotted as arrows showing the correlation between each variable and both PC1 and PC2.
> [!IMPORTANT]
> These arrows only represent the ***direction*** **NOT statistical signficance or strength of the correlation**

<img width="7800" height="5100" alt="ggplot_f3" src="https://github.com/user-attachments/assets/d6d1c8ca-34ff-4d0f-8cd9-3dd6f4d48370" />

## Section 3: Differential Expression Analysis
[DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) is a tool that can complete differential expression analysis (DEA) on the expression matrices to identify differentially expressed genes (DEGs).\
When generating a DESeq2 object the design used by the tool can be altered to include addition variables (i.e. confounding variables) to account for the impact the variable has on the differential expression detected.
```
dds <- DESeq2::DESeqDataSetFromMatrix(countData = <dataset_counts>,
                                      colData = <dataset_samples>,
                                      design = ~ <Confounding_Factor> + condition`
```

The results of DEA can be used to determine how many or which DEGs are overlapping between datasets using packages such as [ggVenn](https://cran.r-project.org/web/packages/ggvenn/index.html).\
Example output from ggVenn:

<img width="430" height="440" alt="Screenshot 2025-08-11 130022" src="https://github.com/user-attachments/assets/3948944a-a1c5-4fd8-b561-efbe7ea12ad0" />


Alternatively, if you have a gene of interest you can plot the normalised read counts of the gene on a box plot across several datasets to look for differences in differential expression between samples and/or datasets (as shown below).

<img width="1212" height="434" alt="Screenshot 2025-08-11 125908" src="https://github.com/user-attachments/assets/fe7736ab-f076-49fb-be10-9ec449ba5f40" />


## Section 4: Functional Annotation
Functional annotation can be done in several ways, but a common one is to visualise enriched [KEGG pathways](http://www.kegg.jp/) using tools like [clusterProfiler](https://www.bioconductor.org/packages//2.10/bioc/html/clusterProfiler.html) and [pathview](https://pathview.r-forge.r-project.org/).\
This enables the visualisation of the pathways that the identified significant DEGs are involved in through plots such as a KEGG dot plot.\
For example:\
This KEGG dot plot shows each enriched KEGG pathway as a dot, where the colour shows how statistically significant the enrichment is and the size relates to how many of those significant DEGs are represented within that pathway. 

<img width="550" height="550" alt="Screenshot 2025-08-11 130539" src="https://github.com/user-attachments/assets/784806c2-dfab-4d26-9bcd-b664fa15e3a9" />


If you find a particular pathway of interest, a KEGG pathway map that shows the up- and down regulation of each gene within that pathway can be generated using the same tools.

Such as the KEGG map below which shows the differential expression of each gene within the [cytokine-cytokine receptor interaction](https://www.kegg.jp/entry/hsa04060) pathway for an individual dataset (down-regulated genes shown in orange and upregulated genes shown in blue)

<img width="1799" height="1161" alt="hsa04060 CTR" src="https://github.com/user-attachments/assets/d9752819-fa3f-47b3-9c5c-0849b54cdce4" />

Alternatively, the genes within a specific pathway can be represented using a mirrored bar plot which separated the genes by direction of change in expression and groups the genes according to signficance of change in expression. 

<img width="1400" height="900" alt="hsa_deg_bar" src="https://github.com/user-attachments/assets/a107c3b0-8576-41cb-b3ef-6298546a112a" />


