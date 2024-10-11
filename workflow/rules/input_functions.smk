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
        #bam_files = {}
        #for progenitor in PROGENITORS:
        #    bam_files[progenitor] = f"results/star/{sample}/{sample}_{progenitor}_aligned.bam"

        #return(bam_files)
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
        output_ref_files[f"{progenitor}"] = f"results/eagle_rc/{sample}/{sample}_classified{index}.ref.bam"

    return output_ref_files


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


def make_eagle_command(type, input, assemblies, params, output):

    command = f"{input['eagle_bin']} --ngi "
    
    if len(sample_files[sample]) == 2:
        command += "--paired "

    for index, progenitor in enumerate(PROGENITORS):
        command += f"--ref{index + 1}={assemblies[progenitor]} --bam{index + 1}={assemblies[progenitor]} " 

    command += f"-o {params['output_prefix']} "

    if type=="WGBS":
        command += "--bs=3 "
    
    command += f"> {output['reads_list']}"

    return command


# Muli QC input
def multiqc_input(type):
    
    input = []

    input.extend(
        expand("results/fastqc/{read_file}_fastqc.zip",read_file=all_read_files)
    )
    
    input.extend(
            expand("results/qualimap/{sample}/{progenitor}", sample=SAMPLES, progenitor=PROGENITORS)     
        )

    input.extend(
            expand("results/eagle_rc/{sample}/{sample}_classified_reads.list", sample=SAMPLES) 
    )

    if type=="WGBS":
        "Something about conversion"

    return input


