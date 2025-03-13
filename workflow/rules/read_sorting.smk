rule install_eagle:
    output:
        eagle_bin=directory("results/eagle_rc/eagle_intallation/"),
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
        "results/logs/{ALIGNER}/renaming/{sample}_{progenitor}.log", 
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -H {input.original_bam} > {output.original_header}
        sed '/^@SQ/ s/SN:/SN:UNIQUENAME_{wildcards.progenitor})_UNIQUENAME_/' {output.original_header} > {output.renamed_header}
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
    run:
        command=make_eagle_command(input, params.assemblies, params, output)
        shell(" && ".join(command))


rule change_sorted_bam_filenames: 
    input:
        log="results/logs/eagle_rc/sorting/{sample}.log",
    log:
        "results/logs/eagle_rc/renaming_files/{sample}.log",
    run:
        command=make_rename_command(wildcards.sample)
        shell(command)


rule restore_chromosome_names_sorted_bams: 
    input:
        #bams=lambda wildcards: get_sorted_bams(wildcards.sample),
        log="results/logs/eagle_rc/renaming_files/{sample}.log",
    log:
        "results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bam_list=$(find results/eagle_rc/{sample}/tmp_renamed/ -name {sample}*.bam)
        for bam in {input.bams}; do
            bad_header=$(echo $bam | 's/.bam/_bad_header.sam/')
            good_header=$(echo $bam | 's/.bam/_good_header.sam/')
            outbam=$(echo $bam | sed 's|/tmp_renamed||')

            samtools view -H $bam > $bad_header
            sed '/^@SQ/ s/SN:UNIQUENAME_*_UNIQUENAME_/SN:/' $bad_header > $good_header
            samtools reheader $good_header $bam > $outbam
            rm $bam $good_header $bad_header
        done

        rm -r results/eagle_rc/{wildcards.sample}/tmp_renamed/
        """