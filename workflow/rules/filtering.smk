rule fastp_se:
    input:
        unpack(lambda wildcards: get_read_files_to_trim(wildcards.sample)),
    output:
        trimmed="results/fastp/{sample}/{sample}_filtered.fastq.gz",
        failed="results/fastp/{sample}{sample}.failed.fastq.gz",
        html="results/fastp/{sample}/{sample}.html",
        json="results/fastp/{sample}/{sample}_se.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: workflow.cores
    params:
        extra=FILTER_PARAMS
    wrapper:
        "v6.2.0/bio/fastp"


rule fastp_pe:
    input:
        unpack(lambda wildcards: get_read_files_to_trim(wildcards.sample)),
    output:
        trimmed=["results/fastp/{sample}/{sample}_R1_filtered.fastq.gz", "results/fastp/{sample}/{sample}_R2_filtered.fastq.gz"],
        unpaired="results/fastp/{sample}/{sample}.unpaired.fastq.gz",
        failed="results/fastp/{sample}/{sample}.failed.fastq.gz",
        html="results/fastp/{sample}/{sample}.html",
        json="results/fastp/{sample}/{sample}_pe.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: workflow.cores
    wrapper:
        "v6.2.0/bio/fastp"