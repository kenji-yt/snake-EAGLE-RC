#### Indexing reference rule ####


# Index genome 1
rule star_index_genome_1:
    input:
        fasta=f"{GENOME_DIR_1}/{GENOME_PARENT_1}.fa",
    output:
        directory(f"{GENOME_DIR_1}/star_genome"),
    message:
        "Testing STAR index"
    threads: workflow.cores
    params:
        extra="",
    log:
        f"logs/star_index_{GENOME_PARENT_1}.log",
    wrapper:
        "v2.3.0/bio/star/index"


# Index genome 2
rule star_index_genome_2:
    input:
        fasta=f"{GENOME_DIR_2}/{GENOME_PARENT_2}.fa",
    output:
        directory(f"{GENOME_DIR_2}/star_genome"),
    message:
        "Testing STAR index"
    threads: workflow.cores
    params:
        extra="",
    log:
        f"logs/star_index_{GENOME_PARENT_1}.log",
    wrapper:
        "v2.3.0/bio/star/index"


#### Alignment rule ####


rule star_pe_multi_1:
    input:
        fq1=f"{RAW_DATA_DIR}/{{sample}}_R1.fastq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=f"{RAW_DATA_DIR}/{{sample}}_R2.fastq.gz",  #optional
        # path to STAR reference genome index
        idx=f"{GENOME_DIR_1}/star_genome",
    output:
        # see STAR manual for additional output files
        aln=f"{OUTPUT_DIR}/star/{{sample}}/1_pe_aligned.bam",
        log=f"{OUTPUT_DIR}/star/{{sample}}/1_Log.out",
        sj=f"{OUTPUT_DIR}/star/{{sample}}/1_SJ.out.tab",
    log:
        "logs/star/{sample}.log",
    params:
        # optional parameters
        extra=f"--outSAMtype BAM SortedByCoordinate",
    threads: workflow.cores
    wrapper:
        "v2.3.0/bio/star/align"


rule star_pe_multi_2:
    input:
        fq1=f"{RAW_DATA_DIR}/{{sample}}_R1.fastq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=f"{RAW_DATA_DIR}/{{sample}}_R2.fastq.gz",  #optional
        # path to STAR reference genome index
        idx=f"{GENOME_DIR_2}/star_genome",
    output:
        # see STAR manual for additional output files
        aln=f"{OUTPUT_DIR}/star/{{sample}}/2_pe_aligned.bam",
        log=f"{OUTPUT_DIR}/star/{{sample}}/2_Log.out",
        sj=f"{OUTPUT_DIR}/star/{{sample}}/2_SJ.out.tab",
    log:
        "logs/star/{sample}.log",
    params:
        # optional parameters
        extra=f"--outSAMtype BAM SortedByCoordinate",
    threads: workflow.cores
    wrapper:
        "v2.3.0/bio/star/align"
