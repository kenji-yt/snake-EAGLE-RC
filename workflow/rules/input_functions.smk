#### Input functions ####

# get fastqc input
def get_fastqc_files(sample):

    if len(sample_files[sample]) == 2:
        return {
            'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}",
            'fq2':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"
        }
    else:
        return {
            'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
        }

# get star input
def get_read_files(sample):

    if len(sample_files[sample]) == 2:
        return {
            'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}",
            'fq2':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"
        }
    else:
        return {
            'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
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
def get_bam_files(type, sample):

    if type=="RNA":
        return {progenitor:f"results/star/{sample}/{sample}_{progenitor}_aligned.bam" for progenitor in PROGENITORS}

    elif type=="DNA": # Unknown yet
        return {progenitor:f"results/??????/{sample}/{sample}_{progenitor}_aligned.bam" for progenitor in PROGENITORS}
    
    elif type=="WGBS": 
        return {progenitor:f"results/bismarks/{sample}/{sample}_{progenitor}_aligned.bam" for progenitor in PROGENITORS}

    else:
        sys.exit(
        f"ERROR: Unknown data type. Name the input directory after the data type: 'DNA', 'RNA' or 'WGBS'."
        )

# get eagle-rc output
def get_eagle_output(sample):
    
    output_ref_files = {}
    for index, progenitor in enumerate(PROGENITORS):
        output_ref_files[progenitor] = f"results/eagle_rc/{sample}/{sample}_classified{index}.ref.bam"

    return output_ref_files

# get assembly
def get_assemblies():
    
    fasta_extensions = ["*.fa", "*.fasta", "*.fna", "*.fq", "fastq"]

    assembly_files = {}
    for progenitor in PROGENITORS:

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


def make_eagle_command(input, params, sample, type):

    command = f"{input.eagle} --ngi "
    
    if len(sample_files[sample]) == 2:
        command += "--paired "

    for index, progenitor in enumerate(PROGENITORS):
        command += f"--ref{index}= {params.{progenitor}} --bam{index}={input.{progenitor}} " 

    command += f"-o {params.output} "

    if type=="WGBS":
        command += "--bs=3 "
    
    command += f"> {params.list}"

    return command


# Muli QC input
def multiqc_input(type):
    
    input = []

    input.extend(
        expand("results/fastqc/{read_file}_fastqc.zip",read_file=all_read_files)
    )
    
    ### ACTUALLY THEY WOULD ALL USE QUALIMAP 
    if type=="RNA":
        input.extend(
            expand("results/qualimap/{sample}/{progenitor}", sample=SAMPLES, progenitor=PROGENITORS)     
        )
    elif type=="DNA":
        "To be determined"
    elif type=="WGBS":
        # I guess :
        input.extend(
            expand("results/bismark/{sample}/{progenitor}", sample=SAMPLES, progenitor=PROGENITORS)     
        )
    else:
        sys.exit(
        f"ERROR: Unknown data type. Name the input directory after the data type: 'DNA', 'RNA' or 'WGBS'."
        )

    return input


