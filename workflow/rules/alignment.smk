rule rename_assemblies:
    input:
        lambda wildcards: get_assembly(wildcards.progenitor),
    output:
        "results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa"
    log:
        "results/logs/assembly_renaming/{progenitor}.log"
    shell:
        """
        sed 's/^>/>UNIQUENAME_{wildcards.progenitor}_UNIQUENAME_/' {input} > {output}
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
    threads: workflow.cores
    wrapper:
        "v4.7.2/bio/bwa-mem2/mem"


rule bwa_mem2_index:
    input:
        "results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa",
    output:
        idx=multiext("results/bwa/{progenitor}/{progenitor}",".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    log:
        "results/logs/bwa/index/{progenitor}.log",
    wrapper:
        "v4.7.2/bio/bwa-mem2/index"


###############################
# WGBS alignment with Bismark #
###############################

rule bismark_alignment:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        genome="results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa",
        bismark_indexes_dir="results/renamed_assemblies/{progenitor}/Bisulfite_Genome",
    output:
        bam="results/bismark/{sample}/{sample}_{progenitor}_aligned.bam",  
        report="results/bismark/{sample}/{sample}_{progenitor}_report.txt",
        
    log:
        "results/logs/bismark/alignment/{sample}_{progenitor}.log"
    threads: 5 # Should be the maximal usage for the alignment. See: https://github.com/FelixKrueger/Bismark/issues/706
    params:
        basename='{sample}_{progenitor}' 
    wrapper:
        "v4.7.2/bio/bismark/bismark"



rule bismark_genome_preparation_fa:
    input:
        "results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa",
    output:
        directory("results/renamed_assemblies/{progenitor}/Bisulfite_Genome")
    log:
        "results/logs/bismark/index/{progenitor}/Bisulfite_Genome.log"
    wrapper:
        "v4.7.2/bio/bismark/bismark_genome_preparation"


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
    message:
        "Aligning RNA seq data with STAR."
    params:
        extra=f"--outSAMtype BAM SortedByCoordinate",
    threads: workflow.cores
    wrapper:
        "v4.7.2/bio/star/align"



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
    wrapper:
        "v4.7.2/bio/star/index"
        