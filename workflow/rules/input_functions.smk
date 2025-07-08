#### Input functions ####

# get filter parameters
def get_filter_params():
    if FILTER_PARAMS in ("No filtering", "Default"):
        return ""
    else:
        return FILTER_PARAMS

# get assembly
def get_assembly(progenitor):
    
    path = os.path.join(f"{INPUT_DIR}/progenitors", progenitor)

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


# get read input for trimming
def get_read_files_to_trim(sample):

    if len(sample_files[sample]) == 2:

        fq1_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
        fq2_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"

        if len(fq1_path) != len(fq2_path):
            error_msg = f"Error: Read filenames for {sample} are of different lengths. Make sure the filenames differ only by one character being either '1' or '2' to indicate read pair entry."
            raise ValueError(error_msg)

        diff_position = None
        for i in range(len(fq1_path)):
            if fq1_path[i] != fq2_path[i]:
                if (fq1_path[i] == '1' and fq2_path[i] == '2') or (fq1_path[i] == '2' and fq2_path[i] == '1'):
                    diff_position = i
                    break
                else:
                    error_msg = f"Error: The only differing character in {sample} read filenames should be '1' or '2'."
                    raise ValueError(error_msg)
        
        return{
            'sample':[fq1_path, fq2_path]
        }
            
    elif len(sample_files[sample]) == 1:
        return {
            'sample':f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
        }
    else:
        error_msg = f"Error: Either more than two or zero read files for sample {sample}."
        raise ValueError(error_msg)



# get read input
def get_read_files(sample):

    if FILTER==True:
        if len(sample_files[sample]) == 2:

            fq1_path=f"results/fastp/{sample}/{sample}_R1_filtered.fastq.gz"
            fq2_path=f"results/fastp/{sample}/{sample}_R2_filtered.fastq.gz"

            if DATA_TYPE=="RNA":
                return {'fq1': fq1_path, 'fq2': fq2_path}

            elif DATA_TYPE=="WGBS":
                return {'fq_1': fq1_path, 'fq_2': fq2_path}

            elif DATA_TYPE=="DNA":
                return {'reads':[fq1_path, fq2_path]}
        
        elif len(sample_files[sample]) == 1:

            fq_path=f"results/fastp/{sample}/{sample}_filtered.fastq.gz"

            if DATA_TYPE=="RNA":
                return {'fq1': fq_path}

            elif DATA_TYPE=="WGBS":
                return {'fq': fq_path}

            elif DATA_TYPE=="DNA":
                return {'reads': fq_path}

        else:
            error_msg = f"Error: Either more than two or zero read files for sample {sample}."
            raise ValueError(error_msg)   

    else:

        if DATA_TYPE=="RNA": 
            if len(sample_files[sample]) == 2:
                fq1_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][0]}"
                fq2_path = f"{INPUT_DIR}/polyploids/{sample}/{sample_files[sample][1]}"

                if len(fq1_path) != len(fq2_path):
                    error_msg = f"Error: Read filenames for {sample} are of different lengths. Make sure the filenames differ only by one character being either '1' or '2' to indicate read pair entry."
                    raise ValueError(error_msg)

                diff_position = None
                for i in range(len(fq1_path)):
                    if fq1_path[i] != fq2_path[i]:
                        if (fq1_path[i] == '1' and fq2_path[i] == '2') or (fq1_path[i] == '2' and fq2_path[i] == '1'):
                            diff_position = i
                            break
                        else:
                            error_msg = f"Error: The only differing character in {sample} read filenames should be '1' or '2'."
                            raise ValueError(error_msg)
                
                if fq1_path[diff_position] == '1':
                    return {'fq1': fq1_path, 'fq2': fq2_path}
                else:
                    return {'fq1': fq2_path, 'fq2': fq1_path}
        
                    
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

                if len(fq1_path) != len(fq2_path):
                    error_msg = f"Error: Read filenames for {sample} are of different lengths. Make sure the filenames differ only by one character being either '1' or '2' to indicate read pair entry."
                    raise ValueError(error_msg)

                diff_position = None
                for i in range(len(fq1_path)):
                    if fq1_path[i] != fq2_path[i]:
                        if (fq1_path[i] == '1' and fq2_path[i] == '2') or (fq1_path[i] == '2' and fq2_path[i] == '1'):
                            diff_position = i
                            break
                        else:
                            error_msg = f"Error: The only differing character in {sample} read filenames should be '1' or '2'."
                            raise ValueError(error_msg)
                
                if fq1_path[diff_position] == '1':
                    return {'fq_1': fq1_path, 'fq_2': fq2_path}
                else:
                    return {'fq_1': fq2_path, 'fq_2': fq1_path} 

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

def check_gtf():

    prog_dir = os.path.join(INPUT_DIR, "progenitors")
    subdirs = [os.path.join(prog_dir, d) for d in os.listdir(prog_dir)
               if os.path.isdir(os.path.join(prog_dir, d))]


    gtf_counts = {}

    for subdir in subdirs:
        gtf_files = [f for f in os.listdir(subdir)
                     if f.endswith(".gtf") and os.path.isfile(os.path.join(subdir, f))]
        gtf_counts[subdir] = len(gtf_files)

    counts = list(gtf_counts.values())

    if all(count == 0 for count in counts):
        return False

    if any(count > 1 for count in counts):
        raise ValueError("One or more progenitor directories contain more than one .gtf file.")

    if any(count == 0 for count in counts):
        raise ValueError("Not every progenitor directory has a .gtf file.")

    return True

# get eagle-rc input
def get_bams(sample):

    return {progenitor:f"results/{ALIGNER}/{sample}/{sample}_{progenitor}_aligned_sorted.bam" for progenitor in PROGENITORS}


# get renamed assemblies dictionary
def get_renamed_assemblies_dict(progenitors): 
    
    assembly_files = {}
    for progenitor in progenitors:

        path = f"results/renamed_assemblies/{progenitor}/renamed_{progenitor}_assembly.fa"
        assembly_files[progenitor] = path

    return assembly_files


def make_eagle_command(input, assemblies, params, output):


    output_directory = f"results/eagle_rc/{params['sample_name']}/tmp_renamed/"

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    if len(PROGENITORS) == 2:

        command = f"{input['eagle_installation']}/eagle-rc --ngi --isc "
        
        if len(sample_files[sample]) == 2:
            command += "--paired "

        if DATA_TYPE=="WGBS":
            command += "--bs=3 "

        if DATA_TYPE=="RNA":
            command += "--splice "


        for index, progenitor in enumerate(PROGENITORS):
            command += f"--ref{index + 1}={assemblies[progenitor]} --bam{index + 1}={input[progenitor]} " 

        command += f"-o {params['output_prefix']} "
    
        
        command += f"> {params['output_prefix']}_reads.list "
        
        check_it=""
        for index, progenitor in enumerate(PROGENITORS):
            check_it += f"md5sum {assemblies[progenitor]}\n" 

        script_content = f"#!/bin/bash\n{check_it}{command}\n"
        
        return script_content

    elif len(PROGENITORS) == 3:

        commands = []
        for i in range(len(PROGENITORS)):
            for j in range(i + 1, len(PROGENITORS)):
                single_command = f"{input['eagle_installation']}/eagle-rc --ngi --isc --listonly \
                    --ref1={assemblies[PROGENITORS[i]]} --bam1={input[PROGENITORS[i]]} \
                    --ref2={assemblies[PROGENITORS[j]]} --bam2={input[PROGENITORS[j]]} "

                if DATA_TYPE=="WGBS":
                    single_command += "--bs=3 "

                if DATA_TYPE=="RNA":
                    single_command += "--splice "

                if len(sample_files[sample]) == 2:
                    single_command += "--paired "

                single_command += f"> {params['output_hexa']}_{PROGENITORS[i]}vs_{PROGENITORS[j]}.list"
                
                commands.append(single_command)


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
            -BD {params['output_hexa']}_{PROGENITORS[1]}vs_{PROGENITORS[2]}.list"
        )

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

        command=" && ".join(commands)

        script_content = f"#!/bin/bash\n{command}\n"

        return script_content

    else:
        error_msg=f"ERROR: Wrong number of progenitors found. Only 2 and 3 allowed. Exiting..."
        raise ValueError(error_msg)


def make_rename_command(sample):

    command=""
    if len(PROGENITORS) == 2:

        for index, progenitor in enumerate(PROGENITORS):
            
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.ref.bam results/eagle_rc/{sample}/tmp_renamed/{sample}_classified_{progenitor}.ref.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming_files/{sample}.log && "
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.mul.bam results/eagle_rc/{sample}/tmp_renamed/{sample}_classified_{progenitor}.mul.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming_files/{sample}.log && "
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.alt.bam results/eagle_rc/{sample}/tmp_renamed/{sample}_classified_{progenitor}.alt.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming_files/{sample}.log && "
            command += f"mv results/eagle_rc/{sample}/tmp_renamed/{sample}_classified{index+1}.unk.bam results/eagle_rc/{sample}/tmp_renamed/{sample}_classified_{progenitor}.unk.bam 2>&1 | tee -a  results/logs/eagle_rc/renaming_files/{sample}.log && "

        command = command.rstrip(' && ')

    elif(PROGENITORS) == 3:

        command="echo Hexaploid: no renaming required."

    return command

def make_qualimap_command(sample, log, threads):
    
    commands = []
    commands.append(f"mkdir -p results/qualimap/{sample} 2> {log}")
    for progenitor in PROGENITORS:

        commands.append(
            f"qualimap bamqc -bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.ref.bam -outdir results/qualimap/{sample}/{progenitor} -outformat PDF:HTML -nt {threads} > {log} 2>&1"
        )

    command = " && ".join(commands)
    return command

def make_RNA_qualimap_command(sample, log):

    commands = []
    commands.append(f"mkdir -p results/qualimap_RNA/{sample} 2> {log}")
    for progenitor in PROGENITORS:

        gtf_file = glob.glob(f"{INPUT_DIR}/progenitors/{progenitor}/*.gtf")[0]

        commands.append(
            f"qualimap rnaseq -bam results/eagle_rc/{sample}/{sample}_classified_{progenitor}.ref.bam -gtf {gtf_file} -outdir results/qualimap_RNA/{sample}/{progenitor} -outformat PDF:HTML > {log} 2>&1"
        )

    command = " && ".join(commands)
    return command

# Muli QC input
def multiqc_input(type):
    
    input = []

    # To define the wildcards
    input.extend(
        expand(f"results/{ALIGNER}" + "/{sample}/{sample}_{progenitor}_aligned_sorted.bam", sample=SAMPLES, progenitor=PROGENITORS) 
    )

    any_gtf = check_gtf()

    if any_gtf == False:
        input.extend(
                expand("results/qualimap/{sample}", sample=SAMPLES)  
            )
    else:
        input.extend(
                    expand("results/qualimap_RNA/{sample}", sample=SAMPLES)     
                )


    if FILTER != False:
        for sample in SAMPLES:
            if len(sample_files[sample]) == 1:
                input.extend([f"results/fastp/{sample}/{sample}_se.fastp.json"])
            elif len(sample_files[sample]) == 2:
                input.extend([f"results/fastp/{sample}/{sample}_pe.fastp.json"])


    return input