rule trimming_se_fastp:
    input:
        unpack(lambda wildcards: get_read_files_to_trim(wildcards.sample)),
    output:
        trimmed="results/trimmed/{sample}/{sample}_trimmed.fastq",
        failed="results/trimmed/{sample}.failed.fastq",
        html="results/report/{sample}.html",
        json="results/report/{sample}.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: workflow.cores
    wrapper:
        "v5.10.0/bio/fastp"


rule trimming_pe_fastp:
    input:
        unpack(lambda wildcards: get_read_files_to_trim(wildcards.sample)),
    output:
        trimmed=["results/trimmed/{sample}/{sample}_R1_trimmed.fastq", "results/trimmed/{sample}/{sample}_R2_trimmed.fastq"],
        unpaired="results/trimmed/{sample}/{sample}.unpaired.fastq",
        failed="results/trimmed/{sample}.failed.fastq",
        html="results/report/{sample}.html",
        json="results/report/{sample}.json"
    log:
        "results/logs/fastp/{sample}.log"
    threads: workflow.cores
    wrapper:
        "v5.10.0/bio/fastp"