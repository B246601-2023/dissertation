#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 定义通用参数
params.output_dir = "results"
params.num_tips = 200
params.sd_min = 0.0
params.sd_max = 1
params.sd_step = 0.1
params.reps = 5
params.seed = 42


process simulateTrees {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}/trees", mode: 'copy'

    input:
    val num_tips 
    val sd_min 
    val sd_max 
    val sd_step 
    val reps 
    val seed 
    val out 

    output:
    path "*" 

    script:
    """
    python3 ${projectDir}/simulate_trees_modify.py --num_tips ${num_tips} --sd_min ${sd_min} --sd_max ${sd_max} --sd_step ${sd_step} --reps ${reps} --out ${out} --seed ${seed}
    """
}

process selectAndCleanTrees {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path tree_files

    output:
    path "selected_trees/*tre" 

    script:
    """
    python3 ${projectDir}/select_trees.py --input ${projectDir}/results/trees --out selected_trees --select trees_selected_trees.txt
    """
}

process rescaleTrees {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path selected_tree_files

    output:
    path "all_trees/*tre" 

    script:
    """
    python3 ${projectDir}/rescale_branch_length.py --input_dir ${projectDir}/results/selected_trees --out all_trees
    """

}

process generateSequences {
    conda '/home/weiwen/envs/tree'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path selected_tree_files

    output:
    path "sequences/*.fasta"

    script:
    """
    mkdir -p sequences
    python3 ${projectDir}/generate_sequences.py --tree ${selected_tree_files} --out "sequences/"
    """
}

// process checkSequences {
//     publishDir "${params.output_dir}/check_plots", mode: 'copy'

//     input:
//     path tree
//     path fasta

//     output:
//     path "*"

//     script:
//     """
//     python3 ${projectDir}/check_seq.py --tree ${tree} --fasta_dir ${projectDir}/results/sequences
//     """
// }

// process extractSequence {
//     conda '/home/weiwen/envs/tree'

//     publishDir "${params.output_dir}", mode: 'copy'

//     input:
//     path fasta_file
    

//     output:
//     path "ref/*.fa"
//     path "sequences_clean/*"

//     script:
//     """
//     python3 ${projectDir}/extract_fasta.py --input ${fasta_file} --out_dir ref --seq_dir sequences_clean --target_sequence root
//     """
// }


workflow {
    num_tips = params.num_tips
    sd_min = params.sd_min
    sd_max = params.sd_max
    sd_step = params.sd_step
    reps = params.reps
    seed = params.seed
    out = "./trees"
    
    tree_files = simulateTrees(
        num_tips=num_tips, 
        sd_min=sd_min, 
        sd_max=sd_max, 
        sd_step=sd_step, 
        reps=reps, 
        seed=seed, 
        out=out
    )
    
    selected_tree_files = selectAndCleanTrees(
        tree_files=tree_files
    )
    
    all_trees = rescaleTrees(
        selected_tree_files = selected_tree_files
    )

    sequences = generateSequences(
        selected_tree_files=selected_tree_files.flatten()
    )

    // checkSequences(
    //     fasta = sequences, tree = selected_tree_files.flatten()
    // )

    // extractSequence(fasta_file = sequences.flatten())

    //align = generateAlignments(fasta_file=sequences.flatten(), ref = refs)

    //vcf = generateVcfs(align = align.flatten())

    //pb = generatePbs(vcf = vcf)
}