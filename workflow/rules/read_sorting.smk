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


rule make_read_sorting_script:
    input:
        unpack(lambda wildcards: get_bams(wildcards.sample)),
        eagle_installation="results/eagle_rc/eagle_intallation",
    output:
        "results/eagle_rc/{sample}/sorting_script.sh"
    params:
        sample_name="{sample}",
        assemblies=get_renamed_assemblies(PROGENITORS),
        output_prefix="results/eagle_rc/{sample}/tmp_renamed/{sample}_classified",
        output_hexa="results/eagle_rc/{sample}/tmp_renamed/{sample}" 
    run:
        script_content=make_eagle_command(input, params.assemblies, params, output)
        script_filename = output[0]
        with open(script_filename, "w") as script_file:
            script_file.write(script_content)


rule read_sorting:
    input:
        script="results/eagle_rc/{sample}/sorting_script.sh"
    log:
        "results/logs/eagle_rc/sorting/{sample}.log",
    conda:
        "../envs/read_sorting.yaml"
    shell:
        "bash {input.script}"



rule change_sorted_bam_filenames_and_delete_renamed_assemblies: 
    input:
        log="results/logs/eagle_rc/sorting/{sample}.log",
    log:
        "results/logs/eagle_rc/renaming_files/{sample}.log",
    run:
        command=make_rename_and_remove_command(wildcards.sample)
        shell(command)


rule restore_chromosome_names_sorted_bams: 
    input:
        log="results/logs/eagle_rc/renaming_files/{sample}.log",
    log:
        "results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bam_list=$(find results/eagle_rc/{wildcards.sample}/tmp_renamed/ -name {wildcards.sample}*.bam)
        for bam in $bam_list; do
            bad_header=$(echo $bam | sed 's/.bam/_bad_header.sam/')
            good_header=$(echo $bam | sed 's/.bam/_good_header.sam/')
            outbam=$(echo $bam | sed 's|/tmp_renamed||')

            samtools view -H $bam > $bad_header
            sed '/^@SQ/ s/SN:UNIQUENAME_*_UNIQUENAME_/SN:/' $bad_header > $good_header
            samtools reheader $good_header $bam > $outbam
            rm $bam $good_header $bad_header
        done

        #rm -r results/eagle_rc/{wildcards.sample}/tmp_renamed/
        """