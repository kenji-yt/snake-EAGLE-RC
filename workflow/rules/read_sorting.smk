
#rule download_eagle:
#    output:
#        eagle_install_dir=directory("/results/eagle_rc/eagle_intallation")
#    log:
#        "results/logs/eagle_rc/download_eagle.log",
#    conda:
#       "envs/build_eagle.yaml"
#    shell:
#        "git clone https://github.com/tony-kuo/eagle.git output.eagle_install_dir"


#rule download_htslib:
#    input:
#        eagle_install_dir="/results/eagle_rc/eagle_intallation"
#    output:
#        htslib_install_dir=directory("/results/eagle_rc/eagle_intallation/htslib"),
#    log:
#       "results/logs/eagle_rc/download_htslib.log",
#    conda:
#        "envs/build_eagle.yaml"
#    shell:
#        "git clone --recursive https://github.com/samtools/htslib.git output.htslib_install_dir"


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
        unpack(lambda wildcards: get_both_bams(ALIGNER, wildcards.sample)),
        eagle_bin="results/eagle_rc/eagle_intallation/eagle-rc",
        #genome1="results/star/{sample}/{sample}_genome1_aligned.bam",
        #genome2="results/star/{sample}/{sample}_genome2_aligned.bam", 
        # FIX get_bam_files function! 
        #reads1=f"{OUTPUT_DIR}/{DATA}_alignment/{{sample}}/1_pe_aligned.bam",
        #reads2=f"{OUTPUT_DIR}/{DATA}_alignment/{{sample}}/2_pe_aligned.bam",
    output:
        reads_list="results/eagle_rc/{sample}/{sample}_classified_reads.list",
        #unpack(lambda wildcards: get_eagle_output(wildcards.sample)),
        #o1=f"results/eagle_rc/{{sample}}/{{sample}}_classified1.ref.bam",
        #o2=f"results/eagle_rc/{{sample}}/{{sample}}_classified2.ref.bam",
    log:
        "results/logs/eagle_rc/sorting/{sample}.log",
    params:
        assemblies=get_assemblies(PROGENITORS),
        #genome1=f"{INPUT_DIR}/progenitors/genome1/g1.fa",
        #genome2=f"{INPUT_DIR}/progenitors/genome2/g2.fa",
        output_prefix="results/eagle_rc/{sample}/{sample}_classified", 
        #output=lambda w, output: os.path.splitext(os.path.splitext(output.o1)[0])[0][
        #    :-1
        #],
    run:
        command=make_eagle_command(DATA_TYPE, input, params.assemblies, params, output)
        shell(command)
        #"{input.eagle_bin} --ngi --paired --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list}"


#rule rename_sorted_reads: 
# SOME RULE TO CHANGE the numbers to actual genome assignment. 
    