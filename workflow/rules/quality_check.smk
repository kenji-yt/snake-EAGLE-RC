#################
## FastQC rule ##
#################

# Quality check for the raw data


rule fastqc:
    input:
        fastq=f"{RAW_DATA_DIR}/{{sample}}_{{extension}}.fastq.gz",
    output:
        html=f"{OUTPUT_DIR}/fastqc/{{sample}}_{{extension}}.html",
        zip=f"{OUTPUT_DIR}/fastqc/{{sample}}_{{extension}}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra="--quiet",
    log:
        f"logs/fastqc/{{sample}}_{{extension}}.log",
    threads: workflow.cores
    resources:
        mem_gb=4,
    wrapper:
        "v2.6.0/bio/fastqc"


###################
## Qualimap rule ##
###################

# Quality check for the alignment with STAR


rule qualimap_1:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam=f"{OUTPUT_DIR}/star/{{sample}}/1_pe_aligned.bam",
    output:
        directory(f"{OUTPUT_DIR}/qualimap/{{sample}}/aligned_1"),
    log:
        f"logs/qualimap/{{sample}}/1_pe_aligned.log",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    params:
        extra=f"-nt {workflow.cores}",
    resources:
        mem_gb=4,
    wrapper:
        "v2.3.2/bio/qualimap/bamqc"


rule qualimap_2:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam=f"{OUTPUT_DIR}/star/{{sample}}/2_pe_aligned.bam",
    output:
        directory(f"{OUTPUT_DIR}/qualimap/{{sample}}/aligned_2"),
    log:
        f"logs/qualimap/{{sample}}/2_pe_aligned.log",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    params:
        extra=f"-nt {workflow.cores}",
    resources:
        mem_gb=4,
    wrapper:
        "v2.3.2/bio/qualimap/bamqc"


###################
#### Multi QC #####
###################

# Multi QC report


rule multiqc_dir:
    input:
        input_files=multiqc_input,
    output:
        out=f"{OUTPUT_DIR}/MultiQC/multiqc_report.html",
    params:
        extra="",  # Optional: extra parameters for multiqc. 
        multiqcdir=lambda w, output: os.path.split(output.out)[0],
    log:
        "logs/multiqc.log",
    wrapper:
        "v2.9.1/bio/multiqc"
