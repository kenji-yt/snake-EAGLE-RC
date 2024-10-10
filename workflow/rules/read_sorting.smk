import os.path

BUILD_PREFIX = ".eagle-build"
INSTALL_PREFIX = ".eagle"

ENV_PATH = os.path.join("../..", "envs", "build_eagle.yaml")

HTSLIB_VERSION = "1.18"
HTSLIB_DIR_NAME = f"htslib-{HTSLIB_VERSION}"
HTSLIB_TAR_NAME = f"{HTSLIB_DIR_NAME}.tar.bz2"
HTSLIB_DIR_PATH = os.path.join(BUILD_PREFIX, HTSLIB_DIR_NAME)
HTSLIB_TAR_PATH = os.path.join(BUILD_PREFIX, HTSLIB_TAR_NAME)
HTSLIB_MK_PATH = os.path.join(HTSLIB_DIR_PATH, "Makefile")

EAGLE = ".eagle/bin/eagle-rc"
EAGLE_VERSION = "1.1.3"
EAGLE_DIR_NAME = f"eagle-{EAGLE_VERSION}"
EAGLE_TAR_NAME = f"{EAGLE_DIR_NAME}.tar.gz"
EAGLE_DIR_PATH = os.path.join(BUILD_PREFIX, EAGLE_DIR_NAME)
EAGLE_TAR_PATH = os.path.join(BUILD_PREFIX, EAGLE_TAR_NAME)
EAGLE_MK_PATH = os.path.join(EAGLE_DIR_PATH, "Makefile")
EAGLE_MK_HTSDIR = os.path.join("..", HTSLIB_DIR_NAME)
EAGLE_MK_PREFIX = os.path.join("..", "..", INSTALL_PREFIX)
EAGLE_BIN_PATH = os.path.join(EAGLE_DIR_PATH, "eagle-rc")
EAGLE_BIN_INSTALL_PATH = os.path.join(INSTALL_PREFIX, "bin", "eagle-rc")


# For whatever reason conda's gcc ignores LIBRARY_PATH
# so we must try and override make's LDFLAGS
# eagle decided to use the non standard name LFLAGS
INCLUDES = os.path.join("$CONDA_PREFIX", "include")
LIBS = os.path.join("$CONDA_PREFIX", "lib")
GCC_PREFIX = "x86_64-conda_cos6-linux-gnu-"
CC = f"{GCC_PREFIX}cc"
AR = f"{GCC_PREFIX}ar"
RANLIB = f"{GCC_PREFIX}ranlib"
HTSLIB_MK_FLAGS = f"LDFLAGS=-L{LIBS}"
EAGLE_MK_FLAGS = f'CC={CC} AR={AR} RANLIB={RANLIB} HTSDIR="{EAGLE_MK_HTSDIR}" PREFIX="{EAGLE_MK_PREFIX}" LFLAGS="-L{EAGLE_MK_HTSDIR} -L{LIBS}" LDLIBS="-lm -lz -llzma -lbz2 -lpthread -lcurl -lssl -lcrypto" MAKE="make {HTSLIB_MK_FLAGS}"'
EAGLE_MK_ENV = f'C_INCLUDE_PATH="{INCLUDES}" LIBRARY_PATH="{LIBS}"'


rule download_eagle:
    output:
        eagle_tar=EAGLE_TAR_PATH,
    params:
        eagle_version=EAGLE_VERSION,
    log:
        "results/logs/eagle_rc/download_eagle.log",
    conda:
        ENV_PATH
    shell:
        "wget https://github.com/tony-kuo/eagle/archive/v{params.eagle_version}.tar.gz -O {output.eagle_tar}"


rule extract_eagle:
    input:
        eagle_tar=EAGLE_TAR_PATH,
    output:
        eagle_mk=EAGLE_MK_PATH,
    params:
        build_prefix=lambda w, input: os.path.split(input.eagle_tar)[0],
    log:
        "results/logs/eagle_rc/extract_eagle.log",
    conda:
        ENV_PATH
    shell:
        "tar -C {params.build_prefix} -vxf {input.eagle_tar}"


rule download_htslib:
    output:
        htslib_tar=HTSLIB_TAR_PATH,
    params:
        htslib_version=HTSLIB_VERSION,
        htslib_tar_name=HTSLIB_TAR_NAME,
    log:
        "results/logs/eagle_rc/download_htslib.log",
    conda:
        ENV_PATH
    shell:
        "wget https://github.com/samtools/htslib/releases/download/{params.htslib_version}/{params.htslib_tar_name} -O {output.htslib_tar}"


rule extract_htslib:
    input:
        htslib_tar=HTSLIB_TAR_PATH,
    output:
        htslib_mk=HTSLIB_MK_PATH,
    params:
        build_prefix=lambda w, input: os.path.split(input.htslib_tar)[0],
    log:
        "results/logs/eagle_rc/extract_htslib.log",
    conda:
        ENV_PATH
    shell:
        "tar -C {params.build_prefix} -vxf {input.htslib_tar}"


rule build_eagle:
    input:
        eagle_mk=EAGLE_MK_PATH,
        htslib_mk=HTSLIB_MK_PATH,
    output:
        eagle_bin=EAGLE_BIN_PATH,
    params:
        eagle_mk_env=EAGLE_MK_ENV,
        eagle_dir_path=lambda w, input: os.path.split(input.eagle_mk)[0],
        eagle_mk_flags=EAGLE_MK_FLAGS,
    log:
        "results/logs/eagle_rc/build_eagle.log",
    conda:
        ENV_PATH
    shell:
        "{params.eagle_mk_env} make -C {params.eagle_dir_path} {params.eagle_mk_flags}"


rule install_eagle:
    input:
        eagle_bin=EAGLE_BIN_PATH,
    output:
        eagle_bin=EAGLE_BIN_INSTALL_PATH,
    params:
        eagle_mk_env=EAGLE_MK_ENV,
        eagle_dir_path=lambda w, input: os.path.split(input.eagle_bin)[0],
        eagle_mk_flags=EAGLE_MK_FLAGS,
    log:
        "results/logs/eagle_rc/install_eagle.log",
    conda:
        ENV_PATH
    shell:
        "{params.eagle_mk_env} make -C {params.eagle_dir_path} install {params.eagle_mk_flags} 2> {log}"


rule read_sorting:
    input:
        eagle_bin=EAGLE,
        unpack(lambda wildcards: get_bam_files(DATA_TYPE, wildcards.sample)), 
        # Make it work for any ploidy!
        #reads1=f"{OUTPUT_DIR}/{DATA}_alignment/{{sample}}/1_pe_aligned.bam",
        #reads2=f"{OUTPUT_DIR}/{DATA}_alignment/{{sample}}/2_pe_aligned.bam",
    output:
        unpack(lambda wildcards: get_eagle_output(wildcards.sample)),
        #o1=f"results/eagle_rc/{{sample}}/{{sample}}_classified1.ref.bam",
        #o2=f"results/eagle_rc/{{sample}}/{{sample}}_classified2.ref.bam",
    log:
        "results/logs/eagle_rc/sorting/{sample}_{progenitor}.log",
    conda:
        ENV_PATH
    params:
        list=f"results/eagle_rc/{sample}/{sample}_classified_reads.list",
        unpack(get_assemblies)
        #genome1=f"{GENOME_DIR_1}/{GENOME_PARENT_1}.fa",
        #genome2=f"{GENOME_DIR_2}/{GENOME_PARENT_2}.fa",
        output=f"results/eagle_rc/{sample}/{sample}_classified", 
        #output=lambda w, output: os.path.splitext(os.path.splitext(output.o1)[0])[0][
        #    :-1
        #],
    shell:
        lambda wildcards: make_eagle_command(input, params, wildcards.sample, DATA_TYPE)
        #"{input.eagle_bin} --ngi --paired --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list}"
