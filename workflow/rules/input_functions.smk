# This file includes all functions to define inputs to other rules


def multiqc_input(wildcards):
    input = []
    if config["RUN_FASTQC"]:
        input.extend(
            expand(
                f"{OUTPUT_DIR}/fastqc/{{sample}}_{{extension}}_fastqc.zip",
                sample=samples.name.values.tolist(),
                extension=["R1", "R2"],
            )
        )
    if config["RUN_STAR"]:
        input.extend(
            expand(
                f"{OUTPUT_DIR}/qualimap/{{sample}}/aligned_{{one_or_two}}",
                sample=samples.name.values.tolist(),
                one_or_two=["1", "2"],
            )
        )
    if config["RUN_EAGLE"]:
        input.extend(
            expand(
                f"{OUTPUT_DIR}/read_sorting/{{sample}}/{{sample}}_classified{{one_or_two}}.ref.bam",
                sample=samples.name.values.tolist(),
                one_or_two=["1", "2"],
            )
        )
    if config["RUN_FEATURECOUNTS"]:
        input.extend(
            expand(
                f"{OUTPUT_DIR}/featureCounts/{{sample}}/{{sample}}_subgenome_{{one_or_two}}{{suffix}}",
                sample=samples.name.values.tolist(),
                one_or_two=["1", "2"],
                suffix=[".featureCounts", ".featureCounts.summary"],
            )
        )
    if config["RUN_EDGER"]:
        input.extend(
            expand(
                f"{OUTPUT_DIR}/edgeR/subgenome_{{one_or_two}}/{{outfile}}",
                one_or_two=["1", "2"],
                outfile=["result_table.txt", "MDS.png", "fc_cpm.png"],
            )
        )

    return input


def edgeR_input(wildcards):
    input = []
    input.extend(
        expand(
            f"{OUTPUT_DIR}/featureCounts/{{sample}}/{{sample}}_subgenome_{{one_or_two}}{{suffix}}",
            sample=samples.name.values.tolist(),
            one_or_two=["1", "2"],
            suffix=[".featureCounts", ".featureCounts.summary"],
        )
    )
    return input
# Special parameters for Feature Count
# def feature_count_params(wildcards):
#    input = []
#    if config["PAIRED_END"]:
#       input.extend("--primary --p")
#      input.extend(config["EXTRA_PARAMS"])
# else:
#    input.extend("--primary")
#   input.extend(config["EXTRA_PARAMS"])
# return " ".join(input)
