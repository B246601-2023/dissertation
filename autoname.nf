#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.output_dir = "results"

process gzip {

    publishDir "${params.output_dir}/pb_gz", mode: 'copy'    

    input:
    path pb_files

    output:
    path "*.gz"

    script:
    """
    gzip -c ${pb_files} > ${pb_files}.gz
    """
}

process autolin {

    conda '/home/weiwen/envs/autolin'

    publishDir "${params.output_dir}/annotated_trees_pb", mode: 'copy'

    input:
    path pb_gz

    output:
    path "*.pb"

    script:
    """
    name=\$(basename ${pb_gz} .pb.gz)
    python3 ${projectDir}/propose_sublineages.py --input ${pb_gz} -o \${name}_annoted.pb  --recursive -t 0 -f 0
    """
    
}

process extract_trees{
    conda '/home/weiwen/envs/usher-env'

    publishDir "${params.output_dir}/annotated_trees_nh", mode: 'copy'

    input:
    path pb_file

    output:
    path "*"

    script:
    """
    name=\$(basename ${pb_file} .pb)
    matUtils extract -i ${pb_file} -t \${name}.nh
    """
}

process extract_txt{
    conda '/home/weiwen/envs/usher-env'

    publishDir "${params.output_dir}/annotated_txt", mode: 'copy'

    input:
    path pb_file

    output:
    path "*.txt"

    script:
    """
    name=\$(basename ${pb_file} .pb)
    matUtils extract -i ${pb_file} -C \${name}.txt
    """
}

process autolin_check{
    conda '/home/weiwen/envs/tree'

    publishDir "${params.output_dir}/autolin_check", mode: 'copy'

    input:
    path txt_file

    output:
    path "*"

    script:
    """
    name=\$(basename ${txt_file} .pb)
    python3 ${projectDir}/autolin_compare.py --input ${txt_file} --output \${name}.csv
    """
}


workflow {
    pb_files = Channel.fromPath("${projectDir}/results/pb/*.pb")
    pb_gz = gzip(pb_files)

    annotated_trees_pb = autolin(pb_gz = pb_gz.flatten())

    nh_trees = extract_trees(pb_file = annotated_trees_pb.flatten())

    annotations = extract_txt(pb_file = annotated_trees_pb.flatten())

    check_results = autolin_check(txt_file = annotations.flatten())
}
