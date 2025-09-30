## Downloading Publicly Available RNA-seq datasets
The following code can be used to download SRA files from the [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/) database.

Navigate to the directory where you want the SRA files downloaded
Ensure you have downloaded the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)  

> [!NOTE]
> The following code can either be submitted using a slurm script or run as an active job in your terminal.

1. Download the *.sra files using the prefetch function of the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)\
`prefetch $accession`

2. Files are downloaded into individual directories so they need to be move out of these subdirectories to current location\
`find . -mindepth 2 -type f -print -exec mv {} . \;`

3. Remove the now empty directories\
`rmdir *`

4. Extract fastq files from the .sra files\
`fastq-dump --split-files *.sra`

5. Delete the sra files so only .fastq files remain\
`rm *.sra`

6. Zip fastq files using [pigz](https://zlib.net/pigz/)\
`pigz *.fastq`

## Generating gene expression matrices using [nf-core/rnaseq](https://nf-co.re/rnaseq/3.19.0/) bioinformatics pipeline
Once the fastq files are processed to a uniform format (e.g. to .fastq.gz files) they can be input into the nf-core/rnaseq pipeline.\
For more information of how this pipeline works, the require inputs and pipeline outputs see their extensive [documentation](https://nf-co.re/rnaseq/3.19.0/docs/usage/)

1. Create a sample sheet/manifest which includes each sample to be run using the following code (assumes paired end reads)
```
find -maxdepth 1 -type f -iname "*_1*" > filenames.txt
cat filenames.txt | cut -d "_" -f1 | cut -d "/" -f2 | sort -n > sample_ID.txt
find "$PWD" -maxdepth 1 -type f -iname "*_1*" | sort -n > forward_reads.txt
find "$PWD" -maxdepth 1 -type f -iname "*_2*" | sort -n > reverse_reads.txt
NROW=$(wc -l < sample_ID.txt)
printf 'auto\n%.0s' $(seq $NROW) > strandedness.txt
paste -d, sample_ID.txt forward_reads.txt reverse_reads.txt strandedness.txt > files_cols.txt
echo -e "sample,fastq_1,fastq_2,strandedness" | cat - files_cols.txt > manifest.csv
rm *.txt
```

2. Load the nf-core module and run the RNA-seq piepline specifying the requirement alignment and quantification tools ([STAR](https://github.com/alexdobin/STAR) [Salmon](https://combine-lab.github.io/salmon/))\
This code can also be run as an active job or submitted as a slurm script
The exact version of the workflow and all included tools can be found in the [nf-core rnaseq software file](nf_core_rnaseq_software_version.yml)
```
module load nextflow/24.04.3 singularity/4.1.0-slurm
nextflow run nf-core/rnaseq --input manifest.csv --outdir rnaseq_results --genome GRCh38 -profile singularity --aligner star_salmon
```
