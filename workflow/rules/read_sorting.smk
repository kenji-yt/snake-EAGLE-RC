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


rule read_sorting:
    input:
        unpack(lambda wildcards: get_both_bams(wildcards.sample)),
        eagle_bin="results/eagle_rc/eagle_intallation/eagle-rc",
    output:
        reads_list="results/eagle_rc/{sample}/{sample}_classified_reads.list",
    log:
        "results/logs/eagle_rc/sorting/{sample}.log",
    params:
        assemblies=get_assemblies(PROGENITORS),
        output_prefix="results/eagle_rc/{sample}/{sample}_classified", 
    run:
        command=make_eagle_command(input, params.assemblies, params, output)
        shell(command)


rule rename_sorted_reads: 
    input:
        reads_list="results/eagle_rc/{sample}/{sample}_classified_reads.list",
    output:
        log="results/logs/eagle_rc/renaming/{sample}.log",
    run:
        command=make_rename_command(wildcards.sample)
        shell(command)
