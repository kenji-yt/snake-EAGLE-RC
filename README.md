# `snake-EAGLE-RC`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


**/!\ NOTE: This tools is only available for Linux /!\\**

A Snakemake workflow around [EAGLE-RC](https://github.com/tony-kuo/eagle?tab=readme-ov-file#eagle-rc), a tool designed to assign sequencing reads to its most likely genome of origin given several possible reference genomes. This is very useful for assigning reads from an allopolyploid to a specific subgenome. 

This workflow (snake-EAGLE-RC) is made to facilitate the read sorting process. It does this by automating the whole process so that you can easily process a large number of samples. The workflow also aims to improve reproducibility by providing a reproducibility report at the end which contains the md5sums of each input and (some) output files and the version (or download date) of each tool used in the analysis. 

*snake-EAGLE-RC* works for RNA, DNA and WGBS data and performs read sorting for up to three subgenomes/progenitors (using --ngi mode of EAGLE-RC). 

## Installation 

- [install Snakemake via Conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
- git clone snake-EAGLE-RC: 
```
git clone https://github.com/kenji-yt/snake-EAGLE-RC.git
```

## Input 

The input to snake-EAGLE-RC is a directory. The name of the directory must be the data type contained inside. Allowed options are: DNA, RNA or WGBS. If your data comes from an rna-seq experiment your input folder must be called 'RNA'. Inside it you should have two directories names 'progenitors' and 'polyploids': 
    - The 'progenitors' directory should contain one subdirectory for each reference genome (each progenitor for allopolyploids). The name of each subdirectory should be the species or genome name and it will appear as such in the output files. Within each subdirectory you should have the corresponding reference genome in fasta/fastq format. Name these files as you wish as long as they have one of the following extensions: "fa","fasta","fq","fna","fastq". In EAGLE-RC, you must ensure that chromosome names are different between subgenomes. **This is not the case in snake-EAGLE-RC.** In other words, you do not need to rename the sequences in your assemblies. 
    - The 'polyploids' directory should contain one subdirectory per sample. The name of each subdirectory should be a unique sample identifier (don't give the same name to different sample directories). Each of these sample directories should contain (appropriately filtered/pre-processed) read data in fasta/fastq format. There should be one file if the data is single-end and two files if paired-end. If paired-end, make sure the filenames are identical expect for a '1' and '2' indicating which side of the pair each file contains. You can input a mix of samples having single-end and paired-end data. 


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


## Analysis 

That's it, you are ready to sort reads using EAGLE-RC. From within the "snake-EAGLE-RC/" directory run:
```
snakemake --use-conda --cores N --config INPUT_DIR='your/input/directory'
```
Make sure to have snakemake make installed, to replace 'your/input/directory' with the path to you input directory and 'N' with the number of cores you wish to allocate to the job. If you installed snakemake in a conda environment make sure to activate it (eg. "conda activate snakemake_env").  

The outputs will now be generated in a results directory within the snake-EAGLE-RC directory. 

#### Quality Check

Quality assessments are performed for the input data ([fastqc](https://github.com/s-andrews/FastQC)) and for the alignments ([qualimap](http://qualimap.conesalab.org/)). You can run the workflow only until these quality assessement steps using "--until fastqc" or "--until qualimap". After verifying the quality, runing the snakemake command will resume the workflow from the quality check step onwards. In the end, all quality checks are reported in 'results/MultiQC/multiqc_report.html'.

#### Trimming

By default the workflow performs trimming using [fastq](https://github.com/OpenGene/fastp). You can set the config argument "TRIMMING=false" to avoid trimming.

#### Additional considerations

Depending on your data you might want to remove low quality reads or checking bisulfite conversion rates independently before running the workflow. snake-EAGLE-RC also does not remove duplicate reads as this is desirable only in some situations. You can judge if you wish to do that once the reads have been sorted or remove them pre-alignment. 
