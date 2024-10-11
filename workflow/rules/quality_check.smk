#################
## FastQC rule ##
#################

rule fastqc:
    input:
        "{read_file}",
    output:
        html="results/fastqc/{read_file}.html",
        zip="results/fastqc/{read_file}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
    log:
        "results/logs/fastqc/{read_file}.log",
    wrapper:
         "v4.0.0/bio/fastqc"


###################
## Qualimap rule ##
###################

rule qualimap_star:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam="results/star/{sample}/{sample}_{progenitor}_aligned.bam",
    output:
        directory("results/qualimap/{sample}/{progenitor}"),
    log:
        "results/logs/qualimap/{sample}/{progenitor}.log",
    #params:
        #extra=f"-nt {workflow.cores}",
    wrapper:
        "v2.3.2/bio/qualimap/bamqc"


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
