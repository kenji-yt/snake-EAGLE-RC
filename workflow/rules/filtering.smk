rule filtering_se_fastp:
    input:
        unpack(lambda wildcards: get_read_files_to_trim(wildcards.sample)),
    output:
        trimmed="results/fastp/{sample}/{sample}_filtered.fastq",
        failed="results/fastp/{sample}{sample}.failed.fastq",
        html="results/fastp/{sample}/{sample}.html",
        json="results/fastp/{sample}/{sample}_se.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: workflow.cores
    params:
        extra=FILTER_PARAMS
    wrapper:
        "v6.2.0/bio/fastp"


rule filtering_pe_fastp:
    input:
        unpack(lambda wildcards: get_read_files_to_trim(wildcards.sample)),
    output:
        trimmed=["results/fastp/{sample}/{sample}_R1_filtered.fastq", "results/fastp/{sample}/{sample}_R2_filtered.fastq"],
        unpaired="results/fastp/{sample}/{sample}.unpaired.fastq",
        failed="results/fastp/{sample}/{sample}.failed.fastq",
        html="results/fastp/{sample}/{sample}.html",
        json="results/fastp/{sample}/{sample}_pe.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: workflow.cores
    wrapper:
        "v6.2.0/bio/fastp"