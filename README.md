## Introduction
This repository contains two versions of a script for completing analyses on recurrent pregnancy loss (RPL) RNA sequencing (RNA-seq) datasets accessed from a public repository, the National Centre for Biotechnology Information (NCBI). 
The first version is the complete script used to analyse and generate the data present in the final paper, while the second version is a simplified script designed with the intention of providing an analysis pipeline for other users completing RNA-seq analysis on RPL datasets. 

## Workflow 
The bioinformatics workflow presented here enables the analysis of RNA sequencing data obtained from human reproductive tissues in RPL research. This pipeline requires a sample sheet containing the sample infromation and gene expression matrices generated using the Salmon tool in the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) bioinformatics pipeline. For more information on how to use the nf-core/rnaseq pipeline including the required imputs and expected outputs, please refer to their [documentation](https://nf-co.re/rnaseq/latest/). 
The processes used to download publicly available high throughput RNA-seq datasets and generate the Salmon gene expression matrices (e.g. counts files) can be found in [Data_Download.md](Data_Preparation.md).

This workflow can be split into different sections to complete the following analyses:
* Section 1: Intialising environment and loading required packages and files 
* Section 2: Principal Component Analysis ([PCAtools](https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html))
  - Section 2.1: Generating PCA objects to be used in Sections 2.2-2.4
  - Section 2.2: Principal Component Retention 
  - Section 2.3: Confounding factor identification using Eigencor plots and Pearson's Correlation coefficients
  - Section 2.4: Generate PCA plots with arrows representing numeric variables ([ggplot2](https://ggplot2.tidyverse.org/))
* Section 3: Differential Expression Analysis ([DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html))
* Section 4: Functional Annotation of [KEGG pathways](http://www.kegg.jp/) ([clusterProfiler](https://www.bioconductor.org/packages//2.10/bioc/html/clusterProfiler.html), [pathview](https://pathview.r-forge.r-project.org/))

## Usage
First, the metadata table must be prepared containing *at least* the information shown below:\
> [!NOTE]
> Make sure the the sample IDs assigned to each sample are __unique__ and __easily identifiable__. If your dataset contains technical replicates (i.e. separates runs from the same sample) they will need to have unique sample IDs (eg. SC_C1r1 and SC_C1r2 etc.)

| SRA Run number | Dataset | Sample ID | Condition |
|:---------------|:-------:|:---------:|----------:|
| SRR####### | SC | SC_C1 | Control |
| SRR####### | SC | SC_C1 | Control |
| SRR####### | SC | SC_RPL1 | RPL |
| SRR####### | SC | SC_RPL2 | RPL |
| SRR####### | RS | RS_C1 | Control |
| SRR####### | RS | RS_C2 | Control |
| SRR####### | RS | RS_RPL1 | RPL |
| SRR####### | RS | RS_RPL2 | RPL |

> [!TIP]
> Add all available metadata to this samples sheet after the columns shown above. This could include information such as patient age, body mass index (BMI) and number of previous miscarriages or live births. _Information such as this can be access using databases like the [NCBI Sequence Read Archive (SRA) Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/)_

## Workflow Outputs
An example of the outputs generated for a set of three datasets using this pipeline can be found in [Outputs.md](Outputs.md). The results generated depend on the information being entered into this workflow and may differ from those generated using other tools due to the differences in analysis methods used by each tool. 


## Acknowledgements
This pipeline was developed during Isabella Brown's Honours Research Project as part of her Bachelor of Advanced Biomedical Science (Honours) course.\
Supervised by Dr Paul Whatmore<sup> 1 </sup> and Dr Kylie Munyard<sup> 1 </sup>\
    <sup> 1 </sup> _Curtin Medical Research Insititute, Curtin Medical School, Curtin Univeristy, Bentley, Australia_




