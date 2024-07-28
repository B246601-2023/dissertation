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

process demix_low{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/demix_results_low", mode: 'copy'
    
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

process demix_high{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/demix_results_high", mode: 'copy'
    
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
    variants_cluster = Channel.fromPath("${projectDir}/results/sample_sets_cluster/*.variants.tsv")
    variants_low = Channel.fromPath("${projectDir}/results/sample_sets_low/*.variants.tsv")
    variants_high = Channel.fromPath("${projectDir}/results/sample_sets_high/*.variants.tsv")

    demix_results_basic = demix_basic(variant = variants_basic)
    demix_results_cluster = demix_cluster(variant = variants_cluster)
    demix_results_low = demix_low(variant = variants_cluster)
    demix_results_high = demix_high(variant = variants_cluster)
}