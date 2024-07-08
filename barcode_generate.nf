#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process barcode_generation {
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/barcodes", mode: 'copy'

    input:
    path vcf

    output:
    path "*.csv", emit: csv
    path "*.feather", emit: feather

    script:
    """
    base_name=\$(basename ${vcf} .vcf)
    pb=/home/weiwen/code/results/annotated_trees_pb/\${base_name}_annoted.pb
    out=\${base_name}_barcode 
    python3 ${projectDir}/convert_barcode_modified.py --input \${pb} --vcf ${vcf} --outname \${out}
    """
}

process sampling{
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/sample_sets", mode: 'copy'

    input:
    path barcode

    output:
    path "*.variants.tsv", emit: variants
    path "*.known_lineages.tsv", emit: known_lineages

    script:
    """
    base_name=\$(basename ${barcode} .csv)
    out_name=\${base_name/_barcode/}
    python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 5 -m 2 -M 5 -f skewed -a 4 -o ./\${out_name}
    """
}

// process delete_empty{
//     publishDir "logs", mode: 'copy'

//     input:
//     path variant

//     output:
//     val signal

//     script:
//     """
//     for $
//     content=\$(cat ${variant})
//     if [ "\${content}" == "empty" ]; then
//         echo "Deleting file because it contains only 'empty'."
//         rm "${variant}"
//     fi
//     """
// }

// process demix{
//     conda '/home/weiwen/envs/usher-env'
//     publishDir "${params.output_dir}/demix_results", mode: 'copy'
    
//     input:
//     path variant
//     path barcode

//     output:
//     path "*.tsv"

//     script:
//     """
//     base_name=\$(basename ${barcode} .feather)
//     depth=\${base_name/_barcode/_depth.tsv}
//     variant=\$(basename ${variant} .variants.tsv)
//     out=\${variant}_demix.tsv
//     if [[ "${variant}" != *empty* ]]; then
//         freyja demix ${variant} ${projectDir}/results/depth/\${depth} --barcodes ${barcode} --output \${out}
//     else
//         touch \${out}
//     fi
//     """
// }

workflow{
    vcf = Channel.fromPath("${projectDir}/results/re_vcf/*.vcf")
    barcodes = barcode_generation(vcf = vcf)
    samples = sampling(barcode = barcodes.csv.flatten())
    // demix_results = demix(variant = samples.variants.flatten(), barcode = barcodes.feather.flatten())
}