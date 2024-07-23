#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 定义通用参数
params.output_dir = "results"

process extractSequence {
    conda '/home/weiwen/envs/tree'

    publishDir "${params.output_dir}", mode: 'copy'
    // errorStrategy 'ignore'
    // when:
    // task.exitStatus == 0

    input:
    path fasta_file
    

    output:
    path "ref/*.fa"
    path "sequences_clean/*", emit: fasta
    path "trees_clean/*"

    script:
    """
    python3 ${projectDir}/extract_fasta.py --input ${fasta_file} --out_dir ref --seq_dir sequences_clean --tree_dir ${projectDir}/results/all_trees
    """
}

process checkSequences {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}/check_plots", mode: 'copy'

    input:
    path tree
    path fasta

    output:
    path "*"

    script:
    """
    base=\$(basename ${fasta} .fa)
    python3 ${projectDir}/check_seq.py --tree ${tree} --fasta_dir ${projectDir}/results/sequences_clean    
    """
}

workflow {
    sequences = Channel.fromPath("${projectDir}/results/sequences/*.fasta")

    sequences_clean = extractSequence(fasta_file = sequences)

    trees = Channel.fromPath("${projectDir}/results/selected_trees/*")
    checkSequences(tree = trees,fasta = sequences_clean.fasta.flatten())

    //align = generateAlignments(fasta_file=sequences.flatten(), ref = refs)

    //vcf = generateVcfs(align = align.flatten())

    //pb = generatePbs(vcf = vcf)
}