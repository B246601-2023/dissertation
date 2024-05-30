#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 定义通用参数
params.output_dir = "results"
params.num_tips = 100
params.sd_min = 0.0
params.sd_max = 1
params.sd_step = 0.1
params.reps = 5
params.seed = 42

process simulateTrees {
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
    python /home/weiwen/code/simulate_trees.py --num_tips ${num_tips} --sd_min ${sd_min} --sd_max ${sd_max} --sd_step ${sd_step} --reps ${reps} --out ${out} --seed ${seed}
    """
}

process selectAndCleanTrees {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path tree_files

    output:
    path "selected_trees/*tre" 

    script:
    """
    python /home/weiwen/code/select_trees.py --input /home/weiwen/code/results/trees --out selected_trees --select trees_selected_trees.txt
    """
}

process rescaleTrees {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path selected_tree_files

    output:
    path "all_trees/*tre" 

    script:
    """
    python /home/weiwen/code/rescale_branch_length.py --input_dir /home/weiwen/code/results/selected_trees --out all_trees
    """

}

process generateSequences {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path selected_tree_files

    output:
    path "sequences/*.fasta"

    script:
    """
    mkdir -p sequences
    for tree in /home/weiwen/code/results/selected_trees/*.tre; do
        python /home/weiwen/code/generate_sequences.py --tree \$tree --out "sequences/"
    done
    """
}

process checkSequences {
    publishDir "${params.output_dir}/check_plots", mode: 'copy'

    input:
    path fasta_files

    output:
    path "*"

    script:
    """
    for tree in /home/weiwen/code/results/selected_trees/*.tre; do
        python /home/weiwen/code/check_seq.py --tree \$tree --fasta_dir /home/weiwen/code/results/sequences
    done
    """
}

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
        selected_tree_files=selected_tree_files
    )

    checkSequences(
        fasta_files = sequences
    )
}
