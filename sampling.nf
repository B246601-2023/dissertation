#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process sampling{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/sample_sets", mode: 'copy'

    input:
    path barcode

    output:
    path "*"

    script:
        """
        base_name=\$(basename ${barcode} .csv)
        out_name=\${base_name/_barcode/}
        python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 5 -m 98 -M 98 -f skewed -a 4 -o ./\${out_name}
        """
}

workflow{
    barcodes = Channel.fromPath("${projectDir}/results/barcodes/*.csv")
    samples = sampling(barcode = barcodes)
}