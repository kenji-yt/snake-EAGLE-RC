# RNA alignment with STAR


rule star_alignment:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        idx = "results/star/{progenitor}",
        
    output:
        # see STAR manual for additional output files
        aln="results/star/{sample}/{progenitor}_aligned.bam",
        #log="results/star/{{sample}}/1_Log.out",#???
        #sj="results/star/{{sample}}/1_SJ.out.tab",#???
    log:
        "results/logs/star/aligment/{sample}_{progenitor}.log",
    message:
        "Aligning RNA seq data with STAR."
    params:
        # optional parameters
        extra=f"--outSAMtype BAM SortedByCoordinate",
    threads: workflow.cores
    wrapper:
        "v4.0.0/bio/star/align"



#### Indexing reference rule ####

rule star_index_genomes:
    input:
        fasta = lambda wildcards: get_assembly(wildcards.progenitor)
    output:
        idx = directory("results/star/{progenitor}"),
    message:
        "Indexing reference genome with STAR."
    threads: workflow.cores
    log:
        "results/logs/star/index/{progenitor}.log",
    wrapper:
        "v4.0.0/bio/star/index"