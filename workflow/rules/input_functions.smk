#### Input functions ####

# get star input
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
def get_input_files(sample):

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
    
    if len(fasta_files) == 0:
        sys.exit(
        f"ERROR: No assembly. No fasta file found in {path}. Exiting..."
        )

    return fasta_files[0]


# Muli QC input
def multiqc_input(type):
    
    input = []

    input.extend(
        expand("results/fastqc/{read_file}_fastqc.zip",read_file=all_read_files)
    )
    
    if type=="RNA":
        input.extend(
            expand("results/qualimap/{sample}/{progenitor}", sample=SAMPLES, progenitor=PROGENITORS)     
        )
    elif type=="DNA":
        "watevaaaa"
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


