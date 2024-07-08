#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process depth_generation{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/depth", mode: 'copy'

    input:
    path fasta_file

    output:
    path "*.tsv"

    script:
    """
    base_name=\$(basename ${fasta_file} .fa)
    out_name=\${base_name/_seqfile/_depth.tsv}
    python3 ${projectDir}/depth_generation.py -fasta ${fasta_file} -depth 1000 -out \${out_name}
    """
}

workflow{
    fasta_file = Channel.fromPath("${projectDir}/results/recombination/*.fa")
    depth_file = depth_generation(fasta_file=fasta_file)
}