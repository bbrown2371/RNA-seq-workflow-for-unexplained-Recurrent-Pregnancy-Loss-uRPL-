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

### Section 2.2: Principal Component (PC) Retention
Generates a named list that contains the retained PCs determined by six different PC retention metrics.\
These can then be used to determine which PC should be used as the limiting PC (i.e. a consensus PC) to identify confounding factors in Section 2.3.\
Depending on the method it may produce a single PC which is the limiting PC or it may return a set of PCs that should be retained according to that metric (see below).
| PC retention metric | Example Output |
| :------- | :------ |
| Horn's Parallel Analysis | Limiting PC (e.g. 3) |
| Scree test (Elbow point) | Limiting PC (e.g. 3) |
| Jolliffe's modification of the Kaiser-Guttman method | Retained PCs (e.g. 1 2 3 4 5) |
| Kaiser-Guttman method | Retained PCs (e.g. 1 2) |
| Cumulative Variance (>80%) | Limiting PC (e.g 8) |
| Modified Broken Stick Method | Retained PCs (e.g. 1 2 16 17 18 19) |

### Section 2.3: Confounding factor identification 
Using the retained PCs determined in Section 2.2, [PCAtools](https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html) can generate eigencor plots to determine which numeric variables are confounding factors.\
For example:\
According to this eigencor plot (assuming the limiting PC is the third PC), the only confounding factor is the number of previous losses ("prev_loss"), which has a significant eigenvalue within PC1.
<br /><br />
<img width="1247" height="864" alt="SC_eigen" src="https://github.com/user-attachments/assets/a8b107be-44d3-4cfa-b577-e4ce3453122e" />
<br /><br />
When an eigencor plot is not possible (i.e. in conditions or groups with fewer than 4 samples), manual Pearson's correlation can be used as follows:
```
# Note: p_list$,dataset_name> represents the PCA object for a specified dataset
pc_df <- p_list$<dataset_name>$rotated
cor.test(pc_df$PC1, p$metadata$age)
```

To check if non-numeric or factor variables are substantially influencing gene expression within a dataset, a PCA plot can be generated where each sample is coloured according to that variable.\
For example:
This plot shows each sample of a dataset coloured according to the reported cycle day (LH+) at endometrial tissue biopsy.
<br /><br />
<img width="2700" height="1800" alt="SC_ggplot_cycle" src="https://github.com/user-attachments/assets/d344913a-3915-4eb3-ad92-3236ebc7be2b" />
<br /><br />
### Section 2.4: Generating PCA plot with arrows for numeric variables 
This process has been designed to allow you to select which numeric clinical or technical variable/s you wish to show on the PCA plot. 

For example:\
These PCA plots show several different clinical variables plotted as arrows indicating the correlation between each variable and both PC1 and PC2.
> [!IMPORTANT]
> These arrows only represent the ***direction*** **NOT statistical signficance or strength** of the correlation
<br /><br />
<img width="7800" height="5100" alt="ggplot_f3" src="https://github.com/user-attachments/assets/0f0f007c-18fb-4d33-b418-48899d2e117f" />
<br /><br />

## Section 3: Differential Expression Analysis

[DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) is a tool that can complete differential expression analysis (DEA) on the expression matrices to identify differentially expressed genes (DEGs).\
When generating a DESeq2 object the design used by the tool can be altered to include addition variables/covariates (i.e. confounding variables) to account for the confounding effect of that variable on differential gene expression.\

The following is an example of the code necessary to run DESeq2 on a single dataset. For a more detailed explanation, see [RNAseq_workflow.R](RNAseq_workflow.R)
```
dds <- DESeq2::DESeqDataSetFromMatrix(countData = <dataset_counts>,
                                      colData = <dataset_samples>,
                                      design = ~ <Confounding_Factor> + condition
```
<br /><br />
The results of DEA can be used to determine how many or which DEGs are overlapping between datasets using packages such as [ggVenn](https://cran.r-project.org/web/packages/ggvenn/index.html).\
Example output from ggVenn:
<br /><br />
<img width="13500" height="3000" alt="ggvenn_f4" src="https://github.com/user-attachments/assets/8704c68c-310b-46a2-aa39-ae0d85725583" />
<br /><br />
Alternatively, if you have a gene of interest you can plot the normalised read counts of the gene on a box plot across several datasets to look for similar differences in differential expression between samples and/or datasets (as shown below).
<br /><br />
<img width="3150" height="3000" alt="CDKN2B" src="https://github.com/user-attachments/assets/0b5fc426-3a5e-4ccb-b501-7f98ad6d3292" />
<br /><br />
## Section 4: Functional Annotation
Functional annotation can be done in several ways, but a common method is to visualise enriched [KEGG pathways](http://www.kegg.jp/) using tools like [clusterProfiler](https://www.bioconductor.org/packages//2.10/bioc/html/clusterProfiler.html) and [pathview](https://pathview.r-forge.r-project.org/).\
This enables the visualisation of the pathways that the significant DEGs are acting within through plots such as a KEGG dot plot.

For example:\
These KEGG dot plots show the top enriched enriched KEGG pathway as dots for each dataset, where the colour shows how statistically significant the enrichment is and the size relates to how many significant DEGs are present within that pathway. 
<br /><br />
<img width="5250" height="3600" alt="kegg_f4" src="https://github.com/user-attachments/assets/a3631ca3-7efa-400b-8df0-fac8429abd54" />
<br /><br />
If you find a particular pathway of interest, a KEGG pathway map that shows the up- and down regulation of each gene within that pathway can be generated using the same tools.

For example:\
The KEGG pathway map below shows the differential expression of each gene within the [cytokine-cytokine receptor interaction](https://www.kegg.jp/entry/hsa04060) pathway for an individual dataset (down-regulated genes shown in orange and upregulated genes shown in blue).\
<br /><br />
<img width="1799" height="1161" alt="hsa04060 CTR" src="https://github.com/user-attachments/assets/d9752819-fa3f-47b3-9c5c-0849b54cdce4" />
<br /><br />

Alternatively, the genes within a specific pathway can be presented using a mirrored bar plot which separates the genes by direction of change in expression (up- or downregulation) and groups the genes according to signficance of the change in expression.\
<br /><br />
<img width="1400" height="900" alt="hsa_deg_bar" src="https://github.com/user-attachments/assets/a107c3b0-8576-41cb-b3ef-6298546a112a" />
<br /><br />

