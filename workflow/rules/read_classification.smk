rule install_eagle:
    output:
        eagle_bin=directory("results/eagle_rc/eagle_installation"),
    log:
        "results/logs/eagle_rc/build_eagle.log",
    params:
        eagle_install_dir="results/eagle_rc/eagle_installation",
        htslib_install_dir="results/eagle_rc/eagle_installation/htslib",
    conda:
        "../envs/build_eagle.yaml"
    shell:
        """
        git clone https://github.com/tony-kuo/eagle.git {params.eagle_install_dir}
        git clone --recursive https://github.com/samtools/htslib.git {params.htslib_install_dir}
        make -C {params.eagle_install_dir}
        """


rule sort_bams:
    input:
        bam=f"results/{ALIGNER}/" + "{sample}/{sample}_{progenitor}_aligned.bam",
    output:
        bam=f"results/{ALIGNER}/" + "{sample}/{sample}_{progenitor}_aligned_sorted.bam",
    log:
        "results/logs/eagle_rc/input_sorting/{sample}/{progenitor}.log",
    conda:
        "../envs/samtools.yaml"
    threads: workflow.cores # adding here so snakemake knows how much this rule uses (idk if necessary)
    shell:
        "samtools sort -T .temp -@ {threads} {input.bam} > {output.bam}"


rule make_read_classification_script:
    input:
        unpack(lambda wildcards: get_bams(wildcards.sample)),
        eagle_installation="results/eagle_rc/eagle_installation",
    output:
        "results/eagle_rc/{sample}/classification_script.sh"
    params:
        sample_name="{sample}",
        assemblies=get_renamed_assemblies_dict(PROGENITORS),
        output_prefix="results/eagle_rc/{sample}/tmp_renamed/{sample}_classified",
        output_hexa="results/eagle_rc/{sample}/tmp_renamed/{sample}",
        classification_log="results/logs/eagle_rc/classification/{sample}.log"
    run:
        script_content=make_eagle_command(input, params.assemblies, params, output)
        script_filename = output[0]
        with open(script_filename, "w") as script_file:
            script_file.write(script_content)


rule read_classification:
    input:
        script="results/eagle_rc/{sample}/classification_script.sh",
    output:
        "results/logs/eagle_rc/classification/{sample}.log",
    params:
        tmp_log="results/logs/eagle_rc/classification/tmp_or_failed_{sample}.log"
    conda:
        "../envs/read_classification.yaml"
    threads: workflow.cores # sad but due to weird https://github.com/tony-kuo/eagle/issues/12#issuecomment-2912831644
    shell:
        """
        bash {input.script} 2>&1 | tee -a {params.tmp_log}
        status=$?
        if [ $status -eq 0 ]; then # Guarantees rerunning following failure. 
            mv {params.tmp_log} {output}
        fi
        exit $status
        """


rule change_classified_bam_filenames: 
    input:
        log="results/logs/eagle_rc/classification/{sample}.log",
    log:
        "results/logs/eagle_rc/renaming_files/{sample}.log",
    run:
        command=make_rename_command(wildcards.sample)
        shell(command)


rule restore_chromosome_names_classified_bams: 
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
            sed '/^@SQ/ s/SN:UNQ_.*_NME_/SN:/' $bad_header > $good_header
            samtools reheader $good_header $bam > $outbam
            rm $bam $good_header $bad_header
        done

        rm -rf results/eagle_rc/{wildcards.sample}/tmp_renamed/

        if [ -d "results/renamed_assemblies" ]; then
            rm -rf "results/renamed_assemblies" 2>&1 | tee -a {log}
        fi
        """
