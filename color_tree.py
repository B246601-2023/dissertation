from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import random
import os
import argparse

def color_subtree(node, color):
    """Colors all subtrees starting from a given node."""
    for n in node.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        nstyle["vt_line_color"] = color
        nstyle["hz_line_color"] = color
        nstyle["vt_line_width"] = 2
        nstyle["hz_line_width"] = 2
        nstyle["vt_line_type"] = 0  # solid
        nstyle["hz_line_type"] = 1  # dashed
        nstyle["fgcolor"] = color
        n.set_style(nstyle)

def apply_color_based_on_labels(tree):
    """Applies colors to nodes based on their labels."""
    for node in tree.traverse():
        if node.name and node.name.startswith("L."):
            color = "#%06x" % random.randint(0, 0xFFFFFF)
            color_subtree(node, color)
            face = TextFace(node.name, fgcolor=color)
            node.add_face(face, column=0, position="branch-right")

def layout(node):
    """Custom layout function to align leaf names to the right."""
    if node.is_leaf():
        color = node.img_style["fgcolor"]
        name_face = TextFace(node.name, fsize=10, fgcolor=color)
        node.add_face(name_face, column=1, position="aligned")

# Set up command line argument parsing
parser = argparse.ArgumentParser(description="Visualize and color a phylogenetic tree.")
parser.add_argument("--input", help="Path to the input tree file")
args = parser.parse_args()

# Load the tree file
tree = Tree(args.input, format=1)

ts = TreeStyle()
apply_color_based_on_labels(tree)

# Visualization settings
ts.min_leaf_separation = 10
ts.draw_guiding_lines = True
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = False
ts.show_scale = False
ts.layout_fn = layout

#tree.show(tree_style=ts)

# Get the base name of the input file for output naming
base_name = os.path.splitext(os.path.basename(args.input))[0]

# Define output file path and format
output_file_path = os.path.join(os.path.dirname(args.input), f"{base_name}.png")
tree.render(output_file_path, w=200, units="mm", tree_style=ts)
print(f"Tree image saved as {output_file_path}")
