# `snake-EAGLE-RC`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow around the eagle-RC tool. The purpose of this workflow is to facilitate the use of [EAGLE-RC](https://github.com/tony-kuo/eagle?tab=readme-ov-file#eagle-rc). 
"EAGLE-RC is a method for classifying whether a read belongs to one genomic hypothesis or another, given a set of genotype differences between them. This can be applicable for determining if reads originate from a specific allele or from a specific homeolog in allopolyploids."



## Usage

The tool offers a minimal amount of filtering post alignment but it is not very flexible in that regard. If you would rather do it yourself and just use the read sorting tool, just create the "results/star" (or/and bismark or XXX depending on data type) directory in the snake-EAGLE-RC directory and add your aligned reads. 

Some rules will automatically determine memory usage but you can set max resources 'snakemake --resources mem_mb=200'.



# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.