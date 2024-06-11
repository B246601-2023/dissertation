#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 定义通用参数
params.output_dir = "results"

process generateAlignments {
    conda '/home/weiwen/envs/usher-env'
    cpus 2
    memory '6 GB'

    publishDir "${params.output_dir}/align", mode: 'copy'

    input:
    path fasta_file 

    output:
    path "*"

    script:
    """
    base_name=\$(basename ${fasta_file} .fa)
    align_fasta="\${base_name}_align.fa"
    root_fasta="\${base_name}_root.fa"
    mafft --thread 2 --auto --keeplength --addfragments ${fasta_file} /home/weiwen/code/results/ref/\${root_fasta} > \${align_fasta}
    """
}

process generateVcfs {
    conda '/home/weiwen/envs/usher-env'

    publishDir "${params.output_dir}/vcf", mode: 'copy'

    input:
    path align  

    output:
    path "*"

    script:
    """
    base_name=\$(basename ${align} .fa)
    vcf_file="\${base_name/_seqfile_align/.vcf}"
    faToVcf ${align} \${vcf_file}
    """
}

process generatePbs {
    conda '/home/weiwen/envs/usher-env'

    publishDir "${params.output_dir}/pb", mode: 'copy'

    input:
    path vcf  

    output:
    path "*"

    script:
    """
    base_name=\$(basename ${vcf} .vcf)
    tree_name="\${base_name}.tre"
    pb_name="\${base_name}.pb"
    usher -t ${projectDir}/results/all_trees/\${tree_name} -v ${vcf} -o \${pb_name}
    """
}

workflow{
    fasta_file = Channel.fromPath("${projectDir}/results/sequences_clean/*")
    align = generateAlignments(fasta_file=fasta_file)

    vcf = generateVcfs(align = align.flatten())

    pb = generatePbs(vcf = vcf)
}