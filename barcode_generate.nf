#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process barcode_generation {
    conda '/home/weiwen/envs/usher-env'
    publishDir "${params.output_dir}/barcodes", mode: 'copy'

    //maxForks 4
    
    input:
    path vcf

    output:
    path "*.csv", emit: csv
    path "*.feather"

    script:
    """
    base_name=\$(basename ${vcf} .vcf)
    base_name="\${base_name//_highre/}"
    pb=${projectDir}/results/annotated_trees/\${base_name}_eddit.pb
    path=${projectDir}/results/autolin_check/\${base_name}_annoted.txt
 
    python3 ${projectDir}/convert_barcode_modified.py --input \${pb} --path \${path} --vcf ${vcf}
    """
}

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
        python3 ${projectDir}/scripts/barcode_sampling.py --barcodes ${barcode} -r 5 -m 2 -M 5 -f skewed -a 4 -o \${out_name}
        """
}

workflow{
    vcf = Channel.fromPath("${projectDir}/results/re_vcf/*_highre.vcf")
    barcodes = barcode_generation(vcf = vcf)
    samples = sampling(barcode = barcodes.csv.flatten())
}