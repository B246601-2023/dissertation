import numpy as np
import dendropy
import toytree
import toyplot.pdf
import toyplot.browser
from dendropy.simulate import treesim
import matplotlib.pyplot as plt


def plot_trees(tree_constant, tree_variable, tree_constant2):
    newick_str_constant = tree_constant.as_string(schema="newick").replace('[&R]', '')
    newick_str_variable = tree_variable.as_string(schema="newick").replace('[&R]', '')
    newick_str_constant2 = tree_constant2.as_string(schema="newick").replace('[&R]', '')


    ttree_constant = toytree.tree(newick_str_constant)
    ttree_variable = toytree.tree(newick_str_variable)
    ttree_constant2 = toytree.tree(newick_str_constant2)
    
    style = {
        "tip_labels_align": True,
       "tip_labels":False,
       "layout": 'd',
        "tip_labels_style": {
            "font-size": "9px"
        },
    }
   
    canvas = toyplot.Canvas(width=1200, height=400)
    axes1 = canvas.cartesian(grid=(1, 3, 0))
    axes2 = canvas.cartesian(grid=(1, 3, 1))
    axes3 = canvas.cartesian(grid=(1, 3, 2))
    
    ttree_constant.draw(axes=axes1, **style)
    ttree_variable.draw(axes=axes2, **style)
    ttree_constant2.draw(axes=axes3, **style)
    
    axes1.label.text = "Constant Rates Tree"
    axes2.label.text = "Variable Rates Tree"

    #toyplot.pdf.render(canvas, out)
    toyplot.browser.show(canvas)

#def simtrees(b):
    # contant rates
tree_constant = treesim.birth_death_tree(
    birth_rate=3,
    death_rate=1,
    birth_rate_sd=0,
    death_rate_sd=0,
    num_extant_tips=50
)

#compare birthrate
tree_constant2 = treesim.birth_death_tree(
    birth_rate=3,
    death_rate=1,
    birth_rate_sd=0,
    death_rate_sd=0,
    num_extant_tips=50
)


    # variable rates
tree_variable = treesim.birth_death_tree(
    birth_rate=3,
    death_rate=1,
    birth_rate_sd=0.5,
    death_rate_sd=0.5,
    num_extant_tips=50
)

#save trees
tree_constant.write(path="tree_constant.tre", schema="newick")
# tree_constant.write(path="tree_constant", schema="newick")


plot_trees(tree_constant, tree_variable, tree_constant2)


