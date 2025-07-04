#########################################
## Let's make a reproducibility report ##
#########################################

input_dir=$1
n_cores=$2
filtering_params=$3
softclipping=$4
workflow_dir=$5
report=results/snake_EAGLE_RC_reproducibility_report.txt
CURRENT_DATETIME=$(date +"%Y-%m-%d %H:%M:%S")

echo "******************" >> "${report}"
echo "* snake-EAGLE-RC *" >> "${report}"
echo "******************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "Reproducibility report for snake-EAGLE-RC." >> "${report}"
echo "Run date & time: ${CURRENT_DATETIME}" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

echo "********************" >> "${report}"
echo "*    Parameters    *" >> "${report}"
echo "********************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "Number of allocated cores: ${n_cores}" >> "${report}"
echo "Input directory: ${input_dir}" >> "${report}"
echo "Filtering parameters: ${filtering_params}" >> "${report}"
echo "Soft clipping: ${softclipping}" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

echo "********************" >> "${report}"
echo "* Operating System *" >> "${report}"
echo "********************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

OS=$(uname -s)

if [ "$OS" == "Linux" ]; then
    # For Linux, try to get version from /etc/os-release
    if [ -f /etc/os-release ]; then
        source /etc/os-release
        echo "Operating System: $NAME" >> "${report}"
        echo "Version: $VERSION" >> "${report}"
    else
        echo "Linux OS (version unknown)" >> "${report}"
    fi
# Assume anything else is macOS
else
    echo "Operating System: macOS"  >> "${report}"
    sw_vers  >> "${report}"
fi

echo "" >> "${report}"
echo "" >> "${report}"
echo "**************" >> "${report}"
echo "* INPUT DATA *" >> "${report}"
echo "**************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

# Loop through each file in the input directory
if [ "$OS" == "Linux" ]; then
    echo "Linux md5sum checksums for the input files" >> "${report}"
    for sub_input_dir in "${input_dir}"/*; do
        if [ -d "${sub_input_dir}" ]; then
            echo $(basename "${sub_input_dir}") >> "${report}"
            for sub_sub in "${sub_input_dir}"/*; do
                if [ -d "${sub_sub}" ]; then
                    echo $(basename "${sub_sub}") >> "${report}"
                    find "${sub_sub}" -type f | xargs -n${n_cores} md5sum | awk '{print $2"\t"$1}' >> "${report}"
                     
                fi
            done
            echo "" >> "${report}"
            echo "" >> "${report}"
        fi
    done
# Assume anything else is macOS
else
    echo "Mac md5 checksums for the input files" >> "${report}"
    for sub_input_dir in "${input_dir}"/*; do
        if [ -d "${sub_input_dir}" ]; then
            echo $(basename "${sub_input_dir}") >> "${report}"
            for sub_sub in "${sub_input_dir}"/*; do
                if [ -d "${sub_sub}" ]; then
                    echo $(basename "${sub_sub}") >> "${report}"
                    find "${sub_sub}" -type f | xargs -n${n_cores} md5 | awk '{print $2"\t"$4}' >> "${report}"
                       
                fi
            done
            echo "" >> "${report}"
            echo "" >> "${report}"
        fi
    done
fi

echo "" >> "${report}"
echo "" >> "${report}"
echo "*********" >> "${report}"
echo "* TOOLS *" >> "${report}"
echo "*********" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"


version_snake_eagle_rc=$(git describe --tags --abbrev=0 | sed 's/v//g')
echo "snake-EAGLE-RC=${version_snake_eagle_rc}" >> "${report}"
eagle_version=$(git -C results/eagle_rc/eagle_installation describe --tags --abbrev=0 | sed 's/v//g')
echo "eagle-rc=${eagle_version}" >> "${report}"

echo "General analysis tools used in their snakemake wrappers:" >> "${report}"
echo "snakemake wrappers release=4.7.2" >> "${report}"
if [ $(basename "${input_dir}")=="RNA" ]; then
    echo "  - star=2.7.11b" >> "${report}"
elif [ $(basename "${input_dir}")=="DNA" ]; then
    echo "  - bwa-mem2=2.2.1" >> "${report}"
elif [ $(basename "${input_dir}")=="WGBS" ]; then
    echo "  - bowtie2=2.5.4" >> "${report}"
    echo "  - bismark=0.24.2" >> "${report}"
fi
echo "snakemake wrappers release=6.2.0" >> "${report}"
echo "  - fastp=0.24.1" >> "${report}"
echo "Tool used for sorting bams:" >> "${report}"
grep 'samtools' workflow/envs/samtools.yaml >> "${report}"
echo "Tool used for quality assessment of classified bam files:" >> "${report}"
grep "qualimap" workflow/envs/qualimap.yaml >> "${report}"
echo "Tools used to install eagle:" >> "${report}"
grep 'git' workflow/envs/build_eagle.yaml >> "${report}"
grep 'make' workflow/envs/build_eagle.yaml >> "${report}"
grep 'zlib' workflow/envs/build_eagle.yaml >> "${report}"
grep 'xz' workflow/envs/build_eagle.yaml >> "${report}"
gcc_version=$(gcc --version | awk '{print $4}' | head -1)
echo "  - gcc='${gcc_version}'" >> "${report}"   
conda list binutils | grep "binutils " | awk '{print "  - "$1"="$2}' >> "${report}" 
echo "Tools used for hexaploid read classification:" >> "${report}"
grep 'python' workflow/envs/read_classification.yaml >> "${report}"
grep 'scipy' workflow/envs/read_classification.yaml >> "${report}"
grep 'numpy' workflow/envs/read_classification.yaml >> "${report}"   

echo "" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "****************" >> "${report}"
echo "* OUTPUT FILES *" >> "${report}"
echo "****************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
if [ "$OS" == "Linux" ]; then
    echo "Linux md5sum checksums for the output files" >> "${report}"
    md5sum results/MultiQC/multiqc_report.html | awk '{print $2"\t"$1}' >> "${report}"  
    echo "" >> "${report}"
    find "results/eagle_rc/" -name "*.ref.bam" | xargs -n${n_cores} md5sum | awk '{print $2"\t"$1}' >> "${report}"
else
    echo "Mac md5 checksums for the input files" >> "${report}"
    md5 results/MultiQC/multiqc_report.html | awk '{print $2"\t"$4}' >> "${report}"  
    echo "" >> "${report}"
    find "results/eagle_rc/" -name "*.ref.bam" | xargs -n${n_cores} md5 | awk '{print $2"\t"$4}' >> "${report}"
fi