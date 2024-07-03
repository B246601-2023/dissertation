#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process recombine {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}/recombination", mode: 'copy'    

    input:
    path tree_file

    output:
    path "*.fa"
    path "*.csv"

    script:
    """
    python3 ${projectDir}/recombination_demo.py --tree ${tree_file} --fasta_dir ${projectDir}/results/sequences_clean/
    """
}

workflow{
    tree_files = Channel.fromPath("${projectDir}/results/autolin_check/modified_tree/*.nh")
    recombination = recombine(tree_file = tree_files)
}
