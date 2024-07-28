#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process sampling_basic{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/sample_sets_basic", mode: 'copy'

    input:
    path barcode

    output:
    path "*"

    script:
        """
        base_name=\$(basename ${barcode} .csv)
        python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 10 -m 50 -M 50 -f skewed -a 4 -o ./\${base_name}
        """
}

process sampling_cluster{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/sample_sets_cluster", mode: 'copy'

    input:
    path barcode

    output:
    path "*"

    script:
        """
        base_name=\$(basename ${barcode} .csv)
        python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 10 -m 50 -M 50 -f skewed -a 4 -c -o ./\${base_name}
        """
}

process sampling_low{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/sample_sets_low", mode: 'copy'

    input:
    path barcode

    output:
    path "*"

    script:
        """
        base_name=\$(basename ${barcode} .csv)
        python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 10 -m 50 -M 50 -f skewed -a 4 -c -o ./\${base_name} -l
        """
}

process sampling_high{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/sample_sets_high", mode: 'copy'

    input:
    path barcode

    output:
    path "*"

    script:
        """
        base_name=\$(basename ${barcode} .csv)
        python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 10 -m 50 -M 50 -f skewed -a 4 -c -o ./\${base_name} -H
        """
}

workflow{
    barcodes = Channel.fromPath("${projectDir}/results/barcodes/*.csv")
    samples_1 = sampling_basic(barcode = barcodes)
    samples_2 = sampling_cluster(barcode = barcodes)
    samples_3 = sampling_low(barcode = barcodes)
    samples_4 = sampling_high(barcode = barcodes)
}