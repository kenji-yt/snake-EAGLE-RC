####################
## Qualimap rules ##
####################

rule sort_ref_bams:
    input:
        classification_log="results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    log:
        sort_ref_log="results/logs/qualimap/{sample}/sorting.log",
    conda:
        "../envs/samtools.yaml"
    threads: workflow.cores
    params:
        threads=workflow.cores,
    run:
        sort_ref_command=make_sort_ref_command(wildcards.sample, log, params.threads)
        shell(sort_ref_command)

        
rule qualimap:
    input:
        sort_ref_log="results/logs/qualimap/{sample}/sorting.log",
    output:
        directory("results/qualimap/{sample}"),
    log:
        "results/logs/qualimap/{sample}/qualimap.log",
    conda:
        "../envs/qualimap.yaml"
    threads: workflow.cores 
    resources:
        mem_mb=1200
    run:
        qualimap_command=make_qualimap_command(wildcards.sample, log, threads, resources.mem_mb)
        shell(qualimap_command)

rule RNA_qualimap:
    input:
        sort_ref_log="results/logs/qualimap/{sample}/sorting.log",
    output:
        directory("results/qualimap_RNA/{sample}"),
    log:
        "results/logs/qualimap_RNA/{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    threads: 1
    run:
        RNA_qualimap_command=make_RNA_qualimap_command(wildcards.sample, log)
        shell(RNA_qualimap_command)
        

###################
#### Multi QC #####
###################

# Multi QC report

rule multiqc_dir:
    input:
        multiqc_input(DATA_TYPE),
    output:
        out = "results/MultiQC/multiqc_report.html", 
    log:
        "results/logs/multiqc.log",
    wrapper:
        "v3.3.3/bio/multiqc"
