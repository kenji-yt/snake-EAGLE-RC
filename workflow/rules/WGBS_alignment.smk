# WGBS alignment with Bismark


rule bismark_alignment:
    input:
        unpack(lambda wildcards: get_read_files(DATA_TYPE, wildcards.sample)),
        genome=lambda wildcards: get_assembly(wildcards.progenitor),
        #bismark_indexes_dir="indexes/{genome}/Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        #unpack(lambda wildcards: get_bismark_out(type, wildcards.sample)), # If I want the option ambiguous and unmapped
        #bam_unmapped_1="bams/{sample}_{genome}_unmapped_reads_1.fq.gz",
        #bam_unmapped_2="bams/{sample}_{genome}_unmapped_reads_2.fq.gz",
        #ambiguous_1="bams/{sample}_{genome}_ambiguous_reads_1.fq.gz",
        #ambiguous_2="bams/{sample}_{genome}_ambiguous_reads_2.fq.gz"
        bam="results/bismark/{sample}/{sample}_{progenitor}_aligned.bam", # NOTE MAYBE THE "pe" or "se" is required in the output definition since working only with basename here. 
        #report="results/bismark/{sample}/{sample}_{progenitor}_report.txt",
        #nucleotide_stats="results/bismark/{sample}/{sample}_{progenitor}.nucleotide_stats.txt",
        
    log:
        "logs/bams/{sample}_{progenitor}.log"
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
        #extra=' --ambiguous --unmapped --nucleotide_coverage',
        basename='{sample}_{progenitor}' # NOTE MAYBE THE "pe" or "se" is required in the output definition since working only with basename here. 
    wrapper:
        "v4.7.1/bio/bismark/bismark"
