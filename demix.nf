#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process demix_basic{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/demix_results_basic", mode: 'copy'
    
    input:
    path variant

    output:
    path "*.tsv"

    script:
    """
    base_name=\$(basename ${variant} .variants.tsv)
    barcode=\${base_name/_rep*/.feather}
    depth=\${barcode/.feather/_depth.tsv}
    out=\${base_name}_demix.tsv
    if [[ "${variant}" != *empty* ]]; then
        freyja demix ${variant} ${projectDir}/results/depth/\${depth} --barcodes ${projectDir}/results/barcodes/\${barcode} --output \${out}
    else
        touch \${out}
    fi
    """
}

process demix_cluster{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/demix_results_cluster", mode: 'copy'
    
    input:
    path variant

    output:
    path "*.tsv"

    script:
    """
    base_name=\$(basename ${variant} .variants.tsv)
    barcode=\${base_name/_rep*/.feather}
    depth=\${barcode/.feather/_depth.tsv}
    out=\${base_name}_demix.tsv
    if [[ "${variant}" != *empty* ]]; then
        freyja demix ${variant} ${projectDir}/results/depth/\${depth} --barcodes ${projectDir}/results/barcodes/\${barcode} --output \${out}
    else
        touch \${out}
    fi
    """
}

workflow{
    variants_basic = Channel.fromPath("${projectDir}/results/sample_sets_basic/*.variants.tsv")
    variants_basic = Channel.fromPath("${projectDir}/results/sample_sets_cluster/*.variants.tsv")

    demix_results_basic = demix(variant = variants_basic)
    demix_results_cluster = demix(variant = variants_cluster)
}