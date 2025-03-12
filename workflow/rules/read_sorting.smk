rule install_eagle:
    output:
        eagle_bin="results/eagle_rc/eagle_intallation/eagle-rc",
    log:
        "results/logs/eagle_rc/build_eagle.log",
    params:
        eagle_install_dir="results/eagle_rc/eagle_intallation",
        htslib_install_dir="results/eagle_rc/eagle_intallation/htslib",
    conda:
        "../envs/build_eagle.yaml"
    shell:
        """
        git clone https://github.com/tony-kuo/eagle.git {params.eagle_install_dir}
        git clone --recursive https://github.com/samtools/htslib.git {params.htslib_install_dir}
        make -C {params.eagle_install_dir}
        """


rule rename_chromosomes:
    input:
        original_bam="results/{ALIGNER}/{sample}/{sample}_{progenitor}_aligned.bam",
    output:
        original_header="results/{ALIGNER}/{sample}/{sample}_{progenitor}_original_header.sam",
        renamed_header="results/{ALIGNER}/{sample}/{sample}_{progenitor}_renamed_header.sam",
        renamed_bam="results/{ALIGNER}/{sample}/{sample}_{progenitor}_renamed.bam",
    log:
        "results/logs/eagle_rc/renaming/{sample}_{progenitor}.log", 
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -H {input.original_bam} > {output.original_header}
        sed '/^@SQ/ s/SN:/SN:{wildcards.progenitor}_/' {output.renamed_header} > {output.renamed_header}
        samtools reheader {output.renamed_header} {input.original_bam} > {output.renamed_bam}
        """


rule read_sorting:
    input:
        unpack(lambda wildcards: get_renamed_bams(wildcards.sample)),
        eagle_installation="results/eagle_rc/eagle_intallation",
    log:
        "results/logs/eagle_rc/sorting/{sample}.log",
    params:
        assemblies=get_assemblies(PROGENITORS),
        output_prefix="results/eagle_rc/{sample}/tmp_renamed/{sample}_classified",
        output_hexa="results/eagle_rc/{sample}/tmp_renamed/{sample}" 
    conda:
        "../envs/hexaploid_sorting.yaml"
    run:
        command=make_eagle_command(input, params.assemblies, params, output)
        shell(" && ".join(command))


rule restore_chromosome_names_sorted_bams: 
    input:
        renamed_chr_sorted_bam=
        log="results/logs/eagle_rc/sorting/{sample}.log",
    output:
        original_chr_name_sorted_bam=
    log:
        "results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    params:
        ref_chromosomes=
        sorted_ref_bam_prefix="results/eagle_rc/{sample}/{sample}",
    conda:
        "../envs/samtools.yaml"
    run:
        """
        samtools view -H {input.original_bam} > {output.original_header}
        sed '/^@SQ/ s/SN:/SN:{wildcards.progenitor}_/' {output.renamed_header} > {output.renamed_header}
        samtools reheader {output.renamed_header} {input.original_bam} > {output.renamed_bam}
        """


rule change_sorted_bam_filenames: 
    input:
        reads_list="results/eagle_rc/{sample}/{sample}_classified_reads.list",
    output:
        log="results/logs/eagle_rc/renaming/{sample}.log",
    run:
        command=make_rename_command(wildcards.sample)
        shell(command)
