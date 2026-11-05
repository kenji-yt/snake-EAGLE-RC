rule rename_assemblies:
    input:
        lambda wildcards: get_assembly(wildcards.progenitor),
    output:
        "results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa"
    log:
        "results/logs/assembly_renaming/{progenitor}.log"
    benchmark:
        "results/benchmarks/rename_assemblies/{progenitor}.tsv",
    shell: 
        """ 
        sed 's/^>/>UNQ_{wildcards.progenitor}_NME_/' {input} > {output}
        """

##########################
# DNA alignment with BWA #
##########################

rule bwa_mem2:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        idx=multiext("results/bwa/{progenitor}/{progenitor}",".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    output:
        "results/bwa/{sample}/{sample}_{progenitor}_aligned.bam",
    log:
        "results/logs/bwa/alignment/{sample}_{progenitor}.log",
    benchmark:
        "results/benchmarks/bwa_mem2/{sample}_{progenitor}.tsv",
    threads: workflow.cores
    wrapper:
        WRAPPER_RELEASE + "/bio/bwa-mem2/mem"


rule bwa_mem2_index:
    input:
        "results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa",
    output:
        idx=multiext("results/bwa/{progenitor}/{progenitor}",".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    log:
        "results/logs/bwa/index/{progenitor}.log",
    benchmark:
        "results/benchmarks/bwa_mem2_index/{progenitor}.tsv",
    wrapper:
        WRAPPER_RELEASE + "/bio/bwa-mem2/index"


###############################
# WGBS alignment with Bismark #
###############################

rule bismark_alignment:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        bismark_indexes_dir="results/renamed_assemblies/{progenitor}/bismark",
    output:
        bam="results/bismark/{sample}/{sample}_{progenitor}_aligned.bam",  
        report="results/bismark/{sample}/{sample}_{progenitor}_report.txt",
        
    log:
        "results/logs/bismark/alignment/{sample}_{progenitor}.log"
    benchmark:
        "results/benchmarks/bismark_alignment/{sample}_{progenitor}.tsv",
    threads: 5 # Should be the maximal usage for the alignment. See: https://github.com/FelixKrueger/Bismark/issues/706
    params:
        extra="--local" if SOFT_CLIP else "",
    wrapper:
        WRAPPER_RELEASE + "/bio/bismark/bismark"



rule bismark_genome_preparation_fa:
    input:
        genome="results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa",
    output:
        bismark_genome_dir=directory("results/renamed_assemblies/{progenitor}/bismark"),
    log:
        "results/logs/bismark/index/{progenitor}/Bisulfite_Genome.log"
    benchmark:
        "results/benchmarks/bismark_genome_preparation_fa/{progenitor}.tsv",
    threads: workflow.cores # the wrapper needs at least 2
    wrapper:
        WRAPPER_RELEASE + "/bio/bismark/bismark_genome_preparation"


###########################
# RNA alignment with STAR #
###########################

rule star_alignment:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        idx = "results/star/{progenitor}/index",
        
    output:
        aln="results/star/{sample}/{sample}_{progenitor}_aligned.bam",
    log:
        "results/logs/star/aligment/{sample}_{progenitor}.log",
    benchmark:
        "results/benchmarks/star_alignment/{sample}_{progenitor}.tsv",
    message:
        "Aligning RNA seq data with STAR."
    params:
        extra=f"--outSAMtype BAM SortedByCoordinate",
    threads: workflow.cores
    wrapper:
        WRAPPER_RELEASE + "/bio/star/align"



rule star_index_genomes:
    input:
        fasta = "results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa"
    output:
        idx = directory("results/star/{progenitor}/index"),
    message:
        "Indexing reference genome with STAR."
    threads: workflow.cores
    log:
        "results/logs/star/index/{progenitor}.log",
    benchmark:
        "results/benchmarks/star_index_genomes/{progenitor}.tsv",
    wrapper:
        WRAPPER_RELEASE + "/bio/star/index"
        
