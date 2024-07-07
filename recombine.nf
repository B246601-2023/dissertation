#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process recombine_generate {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}/recombination", mode: 'copy'    

    input:
    path tree_file

    output:
    path "*.fa", emit: sequences
    path "*.csv"

    script:
    """
    python3 ${projectDir}/recombination_demo.py --tree ${tree_file} --fasta_dir ${projectDir}/results/sequences_clean/
    """
}

process align {
    conda '/home/weiwen/envs/usher-env'
    
    input:
    path fasta_file

    output:
    path "*.fa"

    script:
    """
    name=\$(basename ${fasta_file} .fa)
    outname=\${name/_seqfile/_alignment.fa}
    ${projectDir}/global_align.sh -i ${fasta_file} -o \${outname} -t 4 -r root
    """
}

process convert_tovcf {
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/re_vcf", mode: 'copy'

    input:
    path alignment

    output:
    path "*.vcf"

    script:
    """
    base_name=\$(basename ${alignment} .fa)
    vcf_file="\${base_name/_alignment/.vcf}"
    faToVcf -ref=root ${alignment} \${vcf_file}
    """
}

workflow{
    tree_files = Channel.fromPath("${projectDir}/results/autolin_check/modified_tree/*.nh")
    recombination = recombine_generate(tree_file = tree_files)
    alignments = align(fasta_file = recombination.sequences.flatten())
    vcfs = convert_tovcf(alignment = alignments.flatten())
}
