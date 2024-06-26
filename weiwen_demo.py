import dendropy
import toytree
import toyplot.pdf
from dendropy.simulate import treesim
import matplotlib.pyplot as plt


def plot_trees(tree_constant, tree_variable, out):
    newick_str_constant = tree_constant.as_string(schema="newick").replace('[&R]', '')
    newick_str_variable = tree_variable.as_string(schema="newick").replace('[&R]', '')
    
    ttree_constant = toytree.tree(newick_str_constant)
    ttree_variable = toytree.tree(newick_str_variable)
    
    style = {
        "layout": 'd',
        "tip_labels_align": True,
        "tip_labels_style": {
            "font-size": "9px"
        },
    }
    
    canvas = toyplot.Canvas(width=800, height=400)
    axes1 = canvas.cartesian(grid=(1, 2, 0))
    axes2 = canvas.cartesian(grid=(1, 2, 1))
    
    ttree_constant.draw(axes=axes1, **style)
    ttree_variable.draw(axes=axes2, **style)
    
    axes1.label.text = "Constant Rates Tree"
    axes2.label.text = "Variable Rates Tree"
    
    toyplot.pdf.render(canvas, out)

# contant rates
tree_constant = treesim.birth_death_tree(
    birth_rate=1.0,
    death_rate=0.5,
    birth_rate_sd=0,
    death_rate_sd=0,
    num_extant_tips=50
)

# variable rates
tree_variable = treesim.birth_death_tree(
    birth_rate=1.0,
    death_rate=0.5,
    birth_rate_sd=0.5,
    death_rate_sd=0.5,
    num_extant_tips=50
)

plot_trees(tree_constant, tree_variable, "output_1.pdf")
