#########################################
## Let's make a reproducibility report ##
#########################################

input_dir=$1
n_cores=$2
report=results/snake-EAGLE-RC-reproducibility_report.txt
CURRENT_DATETIME=$(date +"%Y-%m-%d %H:%M:%S")

echo "******************" >> "${report}"
echo "* snake-EAGLE-RC *" >> "${report}"
echo "******************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "Reproducibility report for snake-EAGLE-RC." >> "${report}"
echo "Run date & time: ${CURRENT_DATETIME}" >> "${report}"
echo "Number of allocated cores: ${n_cores}" >> "${report}" 
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
for sub_input_dir in "${input_dir}"/*; do
    if [ -d "${sub_input_dir}" ]; then
        echo $(basename "${sub_input_dir}") >> "${report}"
        for sub_sub in "${sub_input_dir}"/*; do
            if [ -d "${sub_sub}" ]; then
                echo $(basename "${sub_sub}") >> "${report}"
                if [ $(basename "${sub_input_dir}") == "progenitors" ]; then
                    find "${sub_sub}" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.fq" -o -name "*.fastq" \) | \
                    xargs -n${n_cores} md5sum | awk '{print $2"\t"$1}' >> "${report}"   
                else
                    find "${sub_sub}" | xargs -n${n_cores} md5sum | awk '{print $2"\t"$1}' >> "${report}"
                fi   
            fi
        done
    echo "" >> "${report}"
    echo "" >> "${report}"
    fi
done

echo "" >> "${report}"
echo "" >> "${report}"
echo "*********" >> "${report}"
echo "* TOOLS *" >> "${report}"
echo "*********" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"


version_snake_eagle_rc=$(git describe --tags --abbrev=0 | sed 's/v//g')
echo "snake-EAGLE-RC=${version_snake_eagle_rc}" >> "${report}"

eagle_version=$(git -C results/eagle_rc/eagle_intallation describe --tags --abbrev=0 | sed '/v//g')
echo "eagle-rc=${eagle_version}"

echo "fastqc=0.12.1" >> "${report}"
if [ $(basename "${input_dir}")=="RNA" ]; then
    echo "star=2.7.11b" >> "${report}"
elif [ $(basename "${input_dir}")=="DNA" ]; then
    echo "bwa-mem2=2.2.1" >> "${report}"
elif [ $(basename "${input_dir}")=="WGBS" ]; then
    echo "bowtie2=2.5.4" >> "${report}"
    echo "bismark=0.24.2" >> "${report}"
fi
echo $(grep 'samtools' workflow/envs/samtools.yaml | sed 's/- s/s/g') >> "${report}" 
echo $(grep 'git' workflow/envs/build_eagle.yaml | sed 's/- s/s/g') >> "${report}"
echo $(grep 'make' workflow/envs/build_eagle.yaml | sed 's/- s/s/g') >> "${report}"  
echo "qualimap=2.3" >> "${report}"
echo "snakemake wrappers release=4.7.2" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "****************" >> "${report}"
echo "* OUTPUT FILES *" >> "${report}"
echo "****************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
md5sum results/MultiQC/multiqc_report.html | awk '{print $2"\t"$1}' >> "${report}"  
echo "" >> "${report}"
find "results/eagle_rc/" -name "*.ref.bam" | xargs -n${n_cores} md5sum | awk '{print $2"\t"$1}' >> "${report}"
