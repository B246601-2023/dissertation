import os
import time
import random
import argparse
import dendropy
import toytree
import toyplot
import toyplot.pdf
import pandas as pd
import numpy as np
import seaborn as sns
from dendropy.simulate import treesim
import matplotlib.pyplot as plt


def main(params):
    # Parameters
    birth_rate = 1.0
    max_reps_plot = 10
    max_tips_plot = 200
    num_tips = params.num_tips
    sd_range = np.arange(params.sd_min, (params.sd_max + params.sd_step),
                         params.sd_step)
    reps = params.reps
    out = params.out

    # seed random number generator
    if params.seed is not None:
        np.random.seed(params.seed)
        random.seed(params.seed+1)
    else:
        np.random.seed(int(time.time()))
        np.random.seed(int(time.time())+1)

    # make any intermediate output directories needed
    os.makedirs(os.path.dirname(out), exist_ok=True)

    # Simulate trees
    trees = simulate_trees(birth_rate=birth_rate, num_tips=num_tips,
                           reps=reps, sd_range=sd_range, out=out)

    # plot trees
    if reps <= max_reps_plot and num_tips <= max_tips_plot:
        plot_trees(trees, out=out)

    # plot some metrics of tree shape/imbalance
    df=plot_tree_stats(trees, out=out)
    # try to select some trees by colless
    trees = select_trees(df,out)
    # simulate sequences


def simulate_trees(birth_rate=1.0, num_tips=100,
                   sd_range=[0.0], reps=1, out=""):
    """
    Simulates and optionally saves a series of phylogenetic trees with varying
    birth rates (Yule model)

    Parameters:
    - birth_rate (float): Speciation rate.
    - num_tips (int): Number of extant taxa in the tree.
    - sd_range (list of float): Standard deviations for rate variability.
    - reps (int): Number of replicates to simulate for each sd value [def.=1]
    - out (str): Output file prefix [default=""]

    Returns:
    - list of tuples: Each tuple contains a dendropy.Tree object and its
                      associated name.

    Example:
    >>> trees = simulate_trees(birth_rate=1.2,
                               num_extant_tips=50, sd_range=[0.1, 0.2],
                               reps=1, out="sim")
    This will generate and save trees to 'sim_s0.1_r1.tre' and
    'sim_s0.2_r1.tre'.
    """
    trees = []
    for sd in sd_range:
        sd = round(sd, 3)
        for r in range(1, reps + 1):
            valid_tree = False
            while not valid_tree:
                tree = dendropy.simulate.birth_death_tree(
                    birth_rate=birth_rate,
                    death_rate=0.0,
                    birth_rate_sd=sd,
                    num_extant_tips=num_tips,
                    repeat_until_success=True
                )
                # Check for negative branch lengths
                valid_tree = all(edge.length is None or edge.length >= 0 for edge in tree.edges())
                if not valid_tree:  # If any length is negative, repeat simulation
                    continue
            
            name = f"{out}_s{sd}_r{r}"
            filename = f"{name}.tre"
            tree.write(path=filename, schema="newick")
            trees.append((tree, name))
    
    return trees


def plot_tree_stats(trees, out="out"):
    # get tree stats
    data = []
    for tree, name in trees:
        # Calculate tree imbalance metrics
        colless = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)
        sackin = dendropy.calculate.treemeasure.sackin_index(tree)
        treeness = dendropy.calculate.treemeasure.treeness(tree)
        gamma = dendropy.calculate.treemeasure.pybus_harvey_gamma(tree)
        b1 = dendropy.calculate.treemeasure.B1(tree)
        nbar = dendropy.calculate.treemeasure.N_bar(tree)
        parsed_params = parse_name(name)
        data.append({
            'SD': parsed_params['sd'],
            'Replicate': parsed_params['rep'],
            'Colless': colless,
            'Sackin': sackin,
            'Treeness': treeness,
            'Pybus&Harvey Gamma': gamma,
            'B1': b1,
            "Nbar": nbar
        })
    df = pd.DataFrame(data)
    #df.to_csv(f"{out}_statistics.csv")

    # dynamic sized plot 
    metrics = [col for col in df.columns if col not in ['SD', 'Replicate']]
    num_metrics = len(metrics)
    num_cols = 3
    num_rows = (num_metrics + num_cols - 1) // num_cols

    # Create subplots
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols,
                             figsize=(num_cols * 5, num_rows * 5))
    axes = axes.flatten()
    for i, metric in enumerate(metrics):
        sns.boxplot(x='SD', y=metric, hue='SD', data=df, ax=axes[i])
        axes[i].set_title(metric)
        axes[i].set_xlabel('Standard Deviation (SD)')
        axes[i].set_ylabel(metric)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{out}_stat_plots.pdf")
    plt.close()
    return(df)

def select_trees(df, out):
    # sort sd and get the unique values
    sd_values = sorted(df['SD'].unique())
    #print(f"all sd values:\n{sd_values}")
    selected_trees = []
    last_colless = None  # store the colless value of the last selected tree

    for i, sd in enumerate(sd_values):
        # select data for current 'sd' group
        current_group = df[df['SD'] == sd]
        #print(f"current sd group : {current_group}")
        if i == 0:
            # for the first group: select the tree with min colless value
            selected_tree = current_group.loc[current_group['Colless'].idxmin()]
            last_colless = selected_tree['Colless']
        elif i == len(sd_values) - 1:
            # for the last group: select the tree with max colless value
            selected_tree = current_group.loc[current_group['Colless'].idxmax()]
        else:
            # for other groups: select a appropriate one
            next_max = df[df['SD'] == sd_values[i + 1]]['Colless'].max()
            valid_trees = current_group[(current_group['Colless'] > last_colless) & (current_group['Colless'] < next_max)]

            if valid_trees.empty:
                raise ValueError(f"No valid tree found for sd={sd} with Colless greater than {last_colless} and less than {next_max}")

            # select appropriate tree
            selected_tree = valid_trees.sort_values(by='Colless', ascending=True).iloc[0]
            last_colless = selected_tree['Colless']

        # generate tree name and print out
        tree_name = f"{out}_s{(selected_tree['SD'])}_r{int(selected_tree['Replicate'])}.tre"
        selected_trees.append(tree_name)
        print(tree_name)
    return selected_trees




def parse_name(name):
    """
    Utility function to parse parameters from a tree name.

    Parameters:
    - filename (str): The name to parse.

    Returns:
    - dict: A dictionary with keys 'out', 'sd', and 'rep'.
    """
    # split by underscores
    parts = name.replace('__', '_').rstrip('_').split('_')
    # loop through and assign to dict
    params = {}
    for i, part in enumerate(parts):
        if i == 0:
            params['out'] = part
        elif part.startswith('s'):
            params['sd'] = float(part[1:])
        elif part.startswith('r'):
            params['rep'] = int(part[1:])
        else:
            print("Unparsed entry in tree name:", part)
    return params


def plot_trees(trees, out=""):
    """
    Plots and saves a PDF with a grid of cladogram plots

    Parameters:
    - trees (list of tuples): Each tuple contains a DendroPy object and name.
    """
    if not trees:
        return

    # group trees by sd
    trees_by_sd = {}
    for tree, name in trees:
        parsed_params = parse_name(name)
        sd = parsed_params['sd']
        if sd not in trees_by_sd:
            trees_by_sd[sd] = []
        trees_by_sd[sd].append((tree, name))

    # dynamic page size
    canvas_width = 300 * max(len(group) for group in trees_by_sd.values())
    canvas_height = 500 * len(trees_by_sd.keys())

    canvas = toyplot.Canvas(height=canvas_height, width=canvas_width)

    # plot trees in grid
    for i, (sd, grouped_trees) in enumerate(sorted(trees_by_sd.items())):
        for j, (tree, name) in enumerate(grouped_trees):
            axes = canvas.cartesian(grid=(len(trees_by_sd),
                                          len(grouped_trees), i, j))
            newick_str = tree.as_string(schema="newick").replace('[&R]', '')
            ttree = toytree.tree(newick_str)
            ttree.ladderize().draw(axes=axes, tip_labels_align=False,
                                   tip_labels=False)
            axes.label.text = f"SD: {sd}, Replicate {parse_name(name)['rep']}"
    toyplot.pdf.render(canvas, f"{out}_tree_plots.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulate phylogeneties with a Yule variable rate model")
    parser.add_argument('--sd_min', type=float, required=False, default=0.0,
                        help="Minimum standard deviation for birth rate")
    parser.add_argument('--sd_max', type=float, required=False, default=1.0,
                        help="Maximum standard deviation for birth rate")
    parser.add_argument('--sd_step', type=float, required=False, default=0.1,
                        help="Step size for standard deviation range.")
    parser.add_argument('--reps', type=int, required=False, default=10,
                        help="Number of replicates per standard deviation")
    parser.add_argument('--num_tips', type=int, required=False, default=100,
                        help="Number of extant tips in the tree.")
    parser.add_argument('--out', type=str, required=False, default="out",
                        help="Output file prefix.")
    parser.add_argument('--seed', type=int, required=False,
                        help="Seed for random number generator.")
    args = parser.parse_args()
    main(args)
