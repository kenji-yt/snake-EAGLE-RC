# `snake-EAGLE-RC`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


**/!\ NOTE: This tools is only available for Linux /!\\**

A Snakemake workflow around [EAGLE-RC](https://github.com/tony-kuo/eagle?tab=readme-ov-file#eagle-rc), a tool designed to classify sequencing reads to their most likely genome of origin given several possible reference genomes. This is very useful for assigning reads from an allopolyploid to a specific subgenome. 

This workflow (snake-EAGLE-RC) is made to facilitate read classification. It does this by automating the whole pipeline so that you can easily process a large number of samples. The workflow also aims to improve reproducibility by producing a reproducibility report containing the md5sums of each input and (some) output files and the version (or download date) of each tool used in the analysis. 

*snake-EAGLE-RC* works for RNA, DNA and WGBS data and performs read sorting for up to three subgenomes/progenitors using the --ngi mode of EAGLE-RC. 

## Workflow

The workflow performs read classification with EAGLE-RC in "--ngi" (no genotype) mode. This requires alignment of polyploid reads to each subgenome assembly individually. The alignment files along with the assemblies are then fed into EAGLE-RC which classifies assigns each red to a given assembly. Here is a detailed breakdown of the workflow:
- Reads are filtered and trimmed using fastp (optional). A quality report is produced. 
- Renamed assemblies are produced to ensure that the chromosome names differ between each assembly (required by EAGLE-RC).
- Reads are aligned to each renamed assembly. This is done with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) for DNA, [star](https://github.com/alexdobin/STAR) for RNA and [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) for WGBS.  
- The bam files are sorted (required for qualimap). 
- Bam files and renamed assemblies are passed to EAGLE-RC which classified reads in each bam file. See [EAGLE-RC](https://github.com/tony-kuo/eagle?tab=readme-ov-file#eagle-rc) for details on the output produced. 
- Classified bam files are modified to restore the original chromosome names. Renamed assemblies are deleted.
- Qualimap produces a report for the bam files containing reads assigned to the assemblies (".ref"). 


## Installation 

- [install Snakemake via Conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
- git clone snake-EAGLE-RC: 
```
git clone https://github.com/kenji-yt/snake-EAGLE-RC.git
```

## Input 

The input to snake-EAGLE-RC is a directory. The name of the directory must be the data type contained inside. Allowed options are: DNA, RNA or WGBS. So if your data comes from an rna-seq experiment your input folder must be called 'RNA'. 
Inside the input directory you should have two directories names 'progenitors' and 'polyploids': 
- **'progenitors':** This directory should contain one subdirectory for each reference genome (each progenitor for allopolyploids). The name of each subdirectory should be the species or genome name and it will appear as such in the output files. Within each subdirectory you should have the corresponding reference genome in fasta/fastq format. Name these files as you wish as long as they have one of the following extensions: "fa","fasta","fq","fna","fastq". In EAGLE-RC, you must ensure that chromosome names are different between subgenomes. **This is not the case in snake-EAGLE-RC.** In other words, you do not need to rename the sequences in your assemblies. **For RNA** data you can add a ".gtf" file in each progenitor directory. This is not essential for read classification but enables quality check of the output bam using [qualimap](https://github.com/EagleGenomics-cookbooks/QualiMap). In the absence of ".gtf" files, qualimap will still make a report but the results will not be very sensible. **For DNA and WGBS** you must avoid having a ".gtf" file in the progenitor directories as this will run qualimap for RNA-seq instead of qualimap for whole genome sequencing data.    
- **'polyploids':** This directory should contain one subdirectory per sample. The name of each subdirectory should be a unique sample identifier (don't give the same name to different sample directories). Each of these sample directories should only contain sequencing reads in fasta/fastq format (gzipped or not). There should be one file if the data is single-end and two files if paired-end. If paired-end, make sure the filenames are identical expect for a '1' and '2' indicating which side of the pair each file contains. You can input a mix of samples having single-end and paired-end data. 


Your input directory should have the following structure:
```
input_directory/
├── progenitors/
│   ├── species_1
│   │   └── assembly.fa
│   ├── species_2
│   │   └── assembly.fa
│   └── species_3
│       └── assembly.fa
└── polyploids/
    ├── sample_1
    │   ├── reads_pe_1.fastq
    │   └── reads_pe_2.fastq
    ├── sample_2
    │   └── reads_se.fastq
    └── sample_3
        ├── reads_pe_1.fastq
        └── reads_pe_2.fastq
```

As mentionned above, you can add ".gtf" files in the "progenitors/speciesN/" directories if you want appropriate qualimap reports. 

## Usage

You are now ready to classify reads using EAGLE-RC. From within the "snake-EAGLE-RC/" directory run:
```
snakemake --use-conda --cores N --config INPUT_DIR='your/input/directory'
```
Make sure to have snakemake make installed, to replace 'your/input/directory' with the path to your input directory and 'N' with the number of cores you wish to allocate to the job. If you installed snakemake in a conda environment make sure to activate it (eg. "conda activate snakemake_env").  

The outputs will now be generated in a results directory within the snake-EAGLE-RC directory. 

## Quality Check & Filtering

Quality check, filtering and trimming can be performed for the input data using [fastp](https://github.com/OpenGene/fastp). You can run the workflow only until the this step using "--until fastp_se" or "fast_pe" for singe-end and paired-end respectively (Assuming only one data type in your input. If you have a mix of paired-end and single-end this wont work). After filtering and verifying the quality, running the snakemake command will resume the workflow from the fastp step onwards.  You can set the config argument "FILTERING=false" to avoid fastp filtering, trimming and reporting. By default, "FILTERING" is set to "True" such that the workflow performs trimming and filtering with fastp in default mode. If you wish to set your own fastp filtering and trimming parameters simply set the "FILTERING" variable to the flags you wish to pass to fastp (eg. "FILTERING='-q 20 -u 30'"). 

Quality reports for the bam files containing the classified reads are also produced using [qualimap](http://qualimap.conesalab.org/). 

In the end, all quality checks are reported in 'results/MultiQC/multiqc_report.html'.

## Output 

Results will be written to a directory called "results" inside the snake-EAGLE-RC directory. In this directory you will find the following files and directories: 
- eagle_rc: Contains the eagle installation and one directory per sample with the eagle-rc results and a script used to produce these results.
- fastp (if FILTER='True'): Contains one directory per sample with filtered and trimmed read files and quality check reports. 
- qualimap (or qualimap_RNA): Contains one directory per sample containing the output of qualimap for every ".ref." bam files in the 'eagle_rc' directory. 
- star/bismark/bwa: Contains bam files of the reads aligned to each subgenome. 
- logs: Contains logs for each analysis.
- MultiQC: Contains the file "multiqc_report.html" which compiles qualimap and fastp reports. 
- snakemake_EAGLE_RC_reproducibility_report.txt: A text file with details about the input and output files and the tools and parameters used. 

## Additional considerations

Snake-EAGLE-RC also does not remove duplicate reads as this is desirable only in some situations. You can judge if you wish to do that once the reads have been sorted or remove them pre-alignment.

For WGBS data, reads which are soft-clipped are automatically removed. You can tell bismark to leave soft-clipped reads with the config argument "SOFT_CLIP=True". By default, this is set to "SOFT_CLIP=False". 
For bisulfite sequencing data, you also want to check bisulfite conversion rates. This must be done independently of the workflow. 
