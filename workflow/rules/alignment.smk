##########################
# DNA alignment with BWA #
##########################

rule bwa_mem:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        idx=multiext(f"{INPUT_DIR}/progenitors/"+"{progenitor}/{progenitor}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/bwa/{sample}/{sample}_{progenitor}_aligned.bam",
    log:
        "results/logs/bwa/alignment/{sample}_{progenitor}.log",
    params:
        #extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        #sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        #sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        #sort_extra="",  # Extra args for samtools/picard.
    threads: workflow.cores
    wrapper:
        "v4.7.1/bio/bwa/mem"


rule bwa_index:
    input:
        lambda wildcards: get_assembly(wildcards.progenitor),
    output:
        idx=multiext(f"{INPUT_DIR}/progenitors/"+"{progenitor}/{progenitor}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "results/logs/bwa/index/{progenitor}.log",
    wrapper:
        "v4.7.1/bio/bwa/index"


###############################
# WGBS alignment with Bismark #
###############################

rule bismark_alignment:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        genome=lambda wildcards: get_assembly(wildcards.progenitor),
        bismark_indexes_dir=f"{INPUT_DIR}"+"/progenitors/{progenitor}/Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        #unpack(lambda wildcards: get_bismark_out(type, wildcards.sample)), # If I want the option ambiguous and unmapped
        #bam_unmapped_1="bams/{sample}_{genome}_unmapped_reads_1.fq.gz",
        #bam_unmapped_2="bams/{sample}_{genome}_unmapped_reads_2.fq.gz",
        #ambiguous_1="bams/{sample}_{genome}_ambiguous_reads_1.fq.gz",
        #ambiguous_2="bams/{sample}_{genome}_ambiguous_reads_2.fq.gz"
        bam="results/bismark/{sample}/{sample}_{progenitor}_aligned.bam", # NOTE MAYBE THE "pe" or "se" is required in the output definition since working only with basename here. 
        report="results/bismark/{sample}/{sample}_{progenitor}_report.txt",
        #nucleotide_stats="results/bismark/{sample}/{sample}_{progenitor}.nucleotide_stats.txt",
        
    log:
        "results/logs/bismark/alignment/{sample}_{progenitor}.log"
    params:
        # optional params string, e.g: -L32 -N0 -X400 --gzip
        # Useful options to tune:
        # (for bowtie2)
        # -N: The maximum number of mismatches permitted in the "seed", i.e. the first L base pairs
        # of the read (deafault: 1)
        # -L: The "seed length" (deafault: 28)
        # -I: The minimum insert size for valid paired-end alignments. ~ min fragment size filter (for
        # PE reads)
        # -X: The maximum insert size for valid paired-end alignments. ~ max fragment size filter (for
        # PE reads)
        # --gzip: Gzip intermediate fastq files
        # --ambiguous --unmapped
        # -p: bowtie2 parallel execution
        # --multicore: bismark parallel execution
        # --temp_dir: tmp dir for intermediate files instead of output directory
        extra=f'--multicore {workflow.cores}',  #' --ambiguous --unmapped --nucleotide_coverage',
        basename='{sample}_{progenitor}' # NOTE MAYBE THE "pe" or "se" is required in the output definition since working only with basename here. 
    wrapper:
        "v4.7.1/bio/bismark/bismark"



rule bismark_genome_preparation_fa:
    input:
        lambda wildcards: get_assembly(wildcards.progenitor),
    output:
        directory(f"{INPUT_DIR}"+"/progenitors/{progenitor}/Bisulfite_Genome")
    log:
        "results/logs/bismark/index/{progenitor}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        "v4.7.1/bio/bismark/bismark_genome_preparation"


###########################
# RNA alignment with STAR #
###########################

rule star_alignment:
    input:
        unpack(lambda wildcards: get_read_files(wildcards.sample)),
        idx = "results/star/{progenitor}",
        
    output:
        # see STAR manual for additional output files
        aln="results/star/{sample}/{sample}_{progenitor}_aligned.bam",
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
        