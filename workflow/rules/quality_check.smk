####################
## Qualimap rules ##
####################

### There is a RNA seq qualimap...
rule qualimap:
    input:
        bam=f"results/{ALIGNER}/" + "{sample}/{sample}_{progenitor}_aligned_sorted.bam",
    output:
        directory("results/qualimap/{sample}/{progenitor}"),
    log:
        "results/logs/qualimap/{sample}/{progenitor}.log",
    resources:
        mem_mb=lambda wildcard, input: input.size_mb+1000
    threads: workflow.cores 
    params:
        extra=f"-nt {workflow.cores}",
    wrapper:
       "v4.7.2/bio/qualimap/bamqc"


rule sort_bams:
    input:
        bam=f"results/{ALIGNER}/" + "{sample}/{sample}_{progenitor}_aligned.bam",
    output:
        bam=f"results/{ALIGNER}/" + "{sample}/{sample}_{progenitor}_aligned_sorted.bam",
    log:
        "results/logs/qualimap/sorting/{sample}/{progenitor}.log",
    conda:
        "../envs/samtools.yaml"
    threads: workflow.cores # adding here so snakemake knows how much this rule uses (idk if necessary)
    shell:
        "samtools sort -T .temp -@ {threads} {input.bam} > {output.bam}"


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
