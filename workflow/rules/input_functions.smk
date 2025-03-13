#### Input functions ####

# get read input
def get_read_files(sample):

    if DATA_TYPE=="RNA": 
        if len(sample_files[sample]) == 2:
            fq1_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
            fq2_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"
            if "_R1" in sample_files[sample][0]:
                return {'fq1': fq1_path, 'fq2': fq2_path}
            elif "_R1" in sample_files[sample][1]:
                return {'fq1': fq2_path, 'fq2': fq1_path}
            else:
                error_msg = f"Error: Could not determine read pair assignment for sample {sample}. Make sure the read files contain '_R1' and '_R2' in their names."
                raise ValueError(error_msg)
        elif len(sample_files[sample]) == 1:
            return {
                    'fq1':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
            }
        else:
            error_msg = f"Error: Either more than two or zero read files for sample {sample}."
            raise ValueError(error_msg)

    elif DATA_TYPE=="WGBS":
        if len(sample_files[sample]) == 2:
            fq1_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
            fq2_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"
            if "_R1" in sample_files[sample][0]:
                return {'fq_1': fq1_path, 'fq_2': fq2_path}
            elif "_R1" in sample_files[sample][1]:
                return {'fq_1': fq2_path, 'fq_2': fq1_path}
            else:
                error_msg = f"Error: Could not determine read pair assignment for sample {sample}. Make sure the read files contain '_R1' and '_R2' in their names."
                raise ValueError(error_msg)
        elif len(sample_files[sample]) == 1:
            return {
                    'fq':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
            }
        else:
            error_msg = f"Error: Either more than two or zero read files for sample {sample}."
            raise ValueError(error_msg)        

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
        error_msg=f"ERROR: Ambigious assembly. More than one fasta file found in {path}. Exiting..."
        raise ValueError(error_msg)
    
    elif len(fasta_files) == 0:
        error_msg=f"ERROR: No assembly. No fasta file found in {path}. Exiting..."
        raise ValueError(error_msg)

    return fasta_files[0]



# get eagle-rc input
def get_renamed_bams(sample):

    return {progenitor:f"results/{ALIGNER}/{sample}/{sample}_{progenitor}_renamed.bam" for progenitor in PROGENITORS}


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
            error_msg=f"ERROR: Ambigious assembly. More than one fasta file found in {path}. Exiting..."
            raise ValueError(error_msg)
    
        elif len(fasta_files) == 0:
            error_msg=f"ERROR: No assembly. No fasta file found in {path}. Exiting..."
            raise ValueError(error_msg)

        assembly_files[progenitor] = fasta_files[0]

    return assembly_files


def make_eagle_command(input, assemblies, params, output):

    if len(PROGENITORS) == 2:

        command = f"{input['eagle_installation']}/eagle-rc --ngi "
        
        if len(sample_files[sample]) == 2:
            command += "--paired "


        for index, progenitor in enumerate(PROGENITORS):
            command += f"--ref{index + 1}={assemblies[progenitor]} --bam{index + 1}={input[progenitor]} " 

        command += f"-o {params['output_prefix']} "
        
        if DATA_TYPE=="WGBS":
            command += "--bs=3 "

        # Maybe add:
        #if DATA_TYPE=="RNA":
        #    command += "--splice "
        
        command += f"> {params['output_prefix']}_reads.list"

        return [command]

    elif len(PROGENITORS) == 3:

        commands = []
        for i in range(len(PROGENITORS)):
            for j in range(i + 1, len(PROGENITORS)):
                commands.append(
                    f"{input['eagle_installation']}/eagle-rc --ngi --listonly \
                    --ref1={assemblies[PROGENITORS[i]]} --bam1={input[PROGENITORS[i]]} \
                    --ref2={assemblies[PROGENITORS[j]]} --bam2={input[PROGENITORS[j]]} \
                    > {params['output_hexa']}_{PROGENITORS[i]}vs_{PROGENITORS[j]}.list"
                )

        if len(sample_files[sample]) == 2:
            pe_flag = "--pe"
        else:
            pe_flag =""

        # The script calles the three genomes A, B and D; hence this notation here. 
        commands.append(
                    f"python {input['eagle_installation']}/scripts/ref3_ngi_consensus.py \
                    {pe_flag} -u -d -o {params['output_hexa']} \
                    -AB {params['output_hexa']}_{PROGENITORS[0]}vs_{PROGENITORS[1]}.list \
                    -AD {params['output_hexa']}_{PROGENITORS[0]}vs_{PROGENITORS[2]}.list \
                    -BD {params['output_hexa']}_{PROGENITORS[1]}vs_{PROGENITORS[2]}.list")

        commands.append(
            f"{input['eagle_installation']}/eagle-rc \
            --refonly --readlist -a {input[PROGENITORS[0]]} \
            -o {params['output_prefix']}_{PROGENITORS[0]} {params['output_hexa']}.chrA.list"
        )

        commands.append(
            f"{input['eagle_installation']}/eagle-rc \
            --refonly --readlist -a {input[PROGENITORS[1]]} \
            -o {params['output_prefix']}_{PROGENITORS[1]} {params['output_hexa']}.chrB.list"
        )

        commands.append(
            f"{input['eagle_installation']}/eagle-rc \
            --refonly --readlist -a {input[PROGENITORS[2]]} \
            -o {params['output_prefix']}_{PROGENITORS[2]} {params['output_hexa']}.chrD.list"
        )

        return commands

    else:
        error_msg=f"ERROR: Wrong number of progenitors found. Only 2 and 3 allowed. Exiting..."
        raise ValueError(error_msg)



def make_rename_command(sample):

    command=""
    if len(PROGENITORS) == 2:

        for index, progenitor in enumerate(PROGENITORS):
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.ref.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.ref.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.mul.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.mul.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.alt.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.alt.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.unk.bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.unk.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming/{sample}.log && "
    
        command = command.rstrip(' && ')

    elif(PROGENITORS) == 3:

        command="echo Hexaploid: no renaming required."

    return command


def get_sorted_bams(sample):
    
    pattern = f"results/eagle_rc/{sample}/tmp_renamed/{sample}*.bam"  # Adjust the directory if needed

    bam_files = glob.glob(pattern)

    if not bam_files:
        error_msg=f"No BAM files found for sample {sample} in results/eagle_rc/{sample}/tmp_renamed."
        raise ValueError(error_msg)

    return bam_files  



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
            expand("results/logs/eagle_rc/restoring_chr_names/{sample}.log", sample=SAMPLES)     
        )

    return input


