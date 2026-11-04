#########################################
## Let's make a reproducibility report ##
#########################################

input_dir=$1
n_cores=$2
filtering_params=$3
softclipping=$4
envs_dir=$5
pipeline_version=$6
wrapper_release=$7
wrappers=$8
eagle_version=$9
htslib_version=${10}
report=results/snake_EAGLE_RC_reproducibility_report.txt
CURRENT_DATETIME=$(date +"%Y-%m-%d %H:%M:%S")

# Remove renamed assemblies directory if it exists
if [ -d "results/renamed_assemblies" ]; then
    rm -rf "results/renamed_assemblies" 
fi 

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


echo "snake-EAGLE-RC=${pipeline_version}" >> "${report}"
echo "eagle-rc=${eagle_version} (built with htslib=${htslib_version})" >> "${report}"

# Tool versions as pinned in the snakemake-wrappers release used
echo "Tools used through snakemake wrappers (release ${wrapper_release}):" >> "${report}"
for wrapper in ${wrappers}; do
    echo "  ${wrapper}:" >> "${report}"
    env_url="https://raw.githubusercontent.com/snakemake/snakemake-wrappers/${wrapper_release}/bio/${wrapper}/environment.yaml"
    if env_yaml=$(curl -fsSL "${env_url}"); then
        echo "${env_yaml}" | sed -n '/^dependencies:/,$p' | grep -E '^[[:space:]]*- ' | sed -E 's/^[[:space:]]*- */    - /; s/ *= */=/' >> "${report}"
    else
        echo "    (versions not available: could not download ${env_url})" >> "${report}"
    fi
done

# Tool versions as pinned in the workflow's own conda environments
echo "Tools in the workflow's conda environments:" >> "${report}"
for env_file in "${envs_dir}"/*.yaml; do
        purpose=$(sed -n 's/^# Used for: *//p' "${env_file}")
    echo "  ${purpose:-$(basename "${env_file}" .yaml)}:" >> "${report}"
    sed -n '/^dependencies:/,$p' "${env_file}" | grep -E '^[[:space:]]*- ' | sed -E 's/^[[:space:]]*- */    - /; s/[[:space:]]*#.*$//' >> "${report}"
done  

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