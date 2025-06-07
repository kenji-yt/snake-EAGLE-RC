####################
## Qualimap rules ##
####################

rule qualimap:
    input:
        classification_log="results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    output:
        directory("results/qualimap/{sample}"),
    log:
        "results/logs/qualimap/{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    threads: workflow.cores 
    params:
        threads=workflow.cores,
    run:
        qualimap_command=make_qualimap_command(wildcards.sample, log, params.threads)
        shell(qualimap_command)

rule RNA_qualimap:
    input:
        classification_log="results/logs/eagle_rc/restoring_chr_names/{sample}.log",
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
    params:
        multiqcdir = lambda w, output: os.path.split(output.out)[0],
    log:
        "results/logs/multiqc.log",
    wrapper:
        "v2.9.1/bio/multiqc"
