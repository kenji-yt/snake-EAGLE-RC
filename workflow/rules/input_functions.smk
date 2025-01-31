#### Input functions ####

# get read input
def get_read_files(sample):

    if DATA_TYPE=="RNA": 
        if len(sample_files[sample]) == 2:
            return {
                'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}",
                'fq2':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"
            }
        else:
            return {
                'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
            }

    elif DATA_TYPE=="WGBS":
        if len(sample_files[sample]) == 2:
            return {
                'fq_1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}",
                'fq_2':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"
            }
        else:
            return {
                'fq':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
            }

    elif DATA_TYPE=="DNA":
        sample_dir = os.path.join(f"{INPUT_DIR}/polyploids", sample)
        files  = os.listdir(sample_dir)
        file_path = [os.path.join(sample_dir, file) for file in files]
        return {
            'reads':file_path
        }


# get assembly
def get_assembly(progenitor):
    
    path = os.path.join(f"{INPUT_DIR}/progenitors",progenitor)

    fasta_extensions = ["*.fa", "*.fasta", "*.fna", "*.fq", "fastq"]

    fasta_files = []

    for extension in fasta_extensions:
        fasta_files.extend(glob.glob(os.path.join(path, extension)))

    if len(fasta_files) > 1:
        sys.exit(
        f"ERROR: Ambigious assembly. More than one fasta file found in {path}. Exiting..."
        )
    
    elif len(fasta_files) == 0:
        sys.exit(
        f"ERROR: No assembly. No fasta file found in {path}. Exiting..."
        )

    return fasta_files[0]


# get eagle-rc input
def get_both_bams(sample):

    return {progenitor:f"results/{ALIGNER}/{sample}/{sample}_{progenitor}_aligned.bam" for progenitor in PROGENITORS}


# get assembly
def get_assemblies(progenitors): # Here I still pass an argument because of previous error. Also check snakemake tutorial.
    
    fasta_extensions = ["*.fa", "*.fasta", "*.fna", "*.fq", "fastq"]

    assembly_files = {}
    for progenitor in progenitors:

        path = os.path.join(f"{INPUT_DIR}/progenitors",progenitor)

        fasta_files = []

        for extension in fasta_extensions:
            fasta_files.extend(glob.glob(os.path.join(path, extension)))

        if len(fasta_files) > 1:
            sys.exit(
            f"ERROR: Ambigious assembly. More than one fasta file found in {path}. Exiting..."
            )
    
        elif len(fasta_files) == 0:
            sys.exit(
            f"ERROR: No assembly. No fasta file found in {path}. Exiting..."
            )

        assembly_files[progenitor] = fasta_files[0]

    return assembly_files


def make_eagle_command(input, assemblies, params, output):

    command = f"{input['eagle_bin']} --ngi "
    
    if len(sample_files[sample]) == 2:
        command += "--paired "

    for index, progenitor in enumerate(PROGENITORS):
        command += f"--ref{index + 1}={assemblies[progenitor]} --bam{index + 1}={input[progenitor]} " 

    command += f"-o {params['output_prefix']} "
    
    if DATA_TYPE=="WGBS":
        command += "--bs=3 "
    
    command += f"> {output['reads_list']}"

    return command



def make_rename_command(sample):

    with open("results/logs/eagle_rc/renaming/{sample}.log", "w") as file:
        pass

    command=""

    for index, progenitor in enumerate(PROGENITORS):
        command += f"mv results/eagle_rc/{sample}/{sample}_classified{index+1}.ref.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.ref.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
        command += f"mv results/eagle_rc/{sample}/{sample}_classified{index+1}.mul.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.mul.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
        command += f"mv results/eagle_rc/{sample}/{sample}_classified{index+1}.alt.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.alt.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
        command += f"mv results/eagle_rc/{sample}/{sample}_classified{index+1}.unk.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.unk.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
    
    command = command.rstrip(' && ')

    return command



# Muli QC input
def multiqc_input(type):
    
    input = []

    input.extend(
        expand("results/fastqc/{read_file}_fastqc.zip", read_file=all_read_files)
    )
    
    input.extend(
            expand("results/qualimap/{sample}/{progenitor}", sample=SAMPLES, progenitor=PROGENITORS)     
        )
    
    input.extend(
            expand("results/logs/eagle_rc/renaming/{sample}.log", sample=SAMPLES)     
        )

    return input


