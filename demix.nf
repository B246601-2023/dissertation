#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process demix{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/demix_results", mode: 'copy'
    
    input:
    path variant

    output:
    path "*.tsv"

    script:
    """
    base_name=\$(basename ${variant} .variants.tsv)
    barcode=\${base_name/_rep*/_barcode.feather}
    depth=\${barcode/_barcode.feather/_depth.tsv}
    out=\${base_name}_demix.tsv
    if [[ "${variant}" != *empty* ]]; then
        freyja demix ${variant} ${projectDir}/results/depth/\${depth} --barcodes ${projectDir}/results/barcodes/\${barcode} --output \${out}
    else
        touch \${out}
    fi
    """
}

workflow{
    variants = Channel.fromPath("${projectDir}/results/sample_sets/*.variants.tsv")
    demix_results = demix(variant = variants)
}