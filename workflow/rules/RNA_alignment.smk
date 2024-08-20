#### Alignment rule ####


rule star_alignment:
    input:
        ### Make an input function 
        branch(
            expand(file_counts["{sample}"],sample=SAMPLES),
            cases={
                2=expand([fq1=sample_files["{sample}"][0], fq2=sample_files["{sample}"][1]], sample=SAMPLES),
                1=expand([fq1=FQ_DICT["{sample}"][0]], sample=SAMPLES)
            }
        )
        idx = expand("results/star/{progenitor}",progenitor=PROGENITORS),

    output:
        # see STAR manual for additional output files
        aln=expand("results/star/{sample}/{progenitor}_aligned.bam", sample=SAMPLES, progenitor=PROGENITORS),
        #log="results/star/{{sample}}/1_Log.out",#???
        #sj="results/star/{{sample}}/1_SJ.out.tab",#???
    log:
        expand("logs/star/aligment/{sample}.log", sample=SAMPLES),
    message:
        "Aligning RNA seq data with STAR."
    params:
        # optional parameters
        extra=f"--outSAMtype BAM SortedByCoordinate",
    threads: workflow.cores
    wrapper:
        "v2.3.0/bio/star/align"



#### Indexing reference rule ####

rule star_index_genomes:
    input:
        fasta = expand(f"{INPUT_DIR}/progenitors/{{progenitor}}/{{file_name}}.fa",progenitor=PROGENITORS),
    output:
        directory(expand("results/star/{progenitor}",progenitor=PROGENITORS)),
    message:
        "Indexing reference genome with STAR."
    threads: workflow.cores
    log:
        expand("logs/star/index/{progenitor}.log",progenitor=PROGENITORS),
    wrapper:
        "v2.3.0/bio/star/index"
