import os
import re
from io import StringIO
from Bio import Phylo
import dendropy
from dendropy.calculate import treecompare

def extract_true_trees():
    trees = []
    diams = ["0.250000", "0.500000", "1.000000", "2.000000", "3.000000", "4.000000"]
    script_dir = os.path.dirname(os.path.abspath(__file__))

    for diameter in diams:
        diam_trees = []
        for index in range(5):
            tree_filename = f"/TIP_testing_data/taxa50diam{diameter}siteCount800/{index}newick.tree"
            with open(script_dir + tree_filename, "r") as file:

                # each file only contains one tree
                file_tree = file.read().strip()
            diam_trees.append(file_tree)
        trees.append(diam_trees)
    return trees


# Takes in two newick strings and compares the two to see if they represent the same tree, return 1 if they are the
# same and 0 if not. the parameter "branchTrue" is meant to be set if branchhLengths are included in the newick string.
# Default set to false.
def compare_tree(newick1, newick2):

    taxon_namespace = dendropy.TaxonNamespace()

    # first convert the newick strings to get object representations for the trees
    tree1 = dendropy.Tree.get_from_string(newick1, schema="newick", taxon_namespace=taxon_namespace)
    tree2 = dendropy.Tree.get_from_string(newick2, schema="newick", taxon_namespace=taxon_namespace)

    rf_distance = treecompare.unweighted_robinson_foulds_distance(tree1, tree2)

    return rf_distance

# method to find coverage for a given tree file
def calc_coverage(treefile_path, true_tree_newick_string):
    with open(treefile_path, 'r') as tree_file:
        content = tree_file.read()

    newick_pattern = re.compile(r'\([^;]*;\s*')
    trees = newick_pattern.findall(content)

    count = 0.0
    for tree in trees:
        count += compare_tree(norm_newick(tree), norm_newick(true_tree_newick_string))
    return count / len(trees)

# helper method to convert Taxa_number_x to Taxax to simplify things. rb naming convention
# can be slightly different
def norm_newick(newick_string): 
    new_string = re.sub(r'Taxa_number_(\d+)', r'Taxa\1', newick_string)
    new_string = re.sub(r':\d*\.?\d*([eE][-+]?\d+)?', ':0.0', new_string)
    return re.sub(r'\[.*?\]', '', new_string)

# run a script to calculate coverage
trees = extract_true_trees()
tree_per_diam = len(trees[0])
diams = ["0.250000", "0.500000", "1.000000", "2.000000", "3.000000", "4.000000"]
script_dir = os.path.dirname(os.path.abspath(__file__))

# lets iterate through the different tree diameters
for index, diameter in enumerate(diams):
    diam_dir = f"taxa50diam{diameter}siteCount800"
    new_model_diam_dir = "/output" + "/" + diam_dir
    rb_model_diam_dir = "/rb_output" + "/" + diam_dir

    # now iterate through the things to check
    for tree_counter in range(tree_per_diam):
        tip_cur_tree_file = f"/TIP{tree_counter}tipSequence.tree"
        rb_cur_tree_file = f"/rb{tree_counter}tipSequence.tree"
        # print(trees[index][tree_counter])
        # print(script_dir + rb_model_diam_dir + rb_cur_tree_file)
        # print(script_dir + new_model_diam_dir + tip_cur_tree_file)
        TIP_coverage = calc_coverage(script_dir + new_model_diam_dir + tip_cur_tree_file, trees[index][tree_counter])
        rb_coverage = calc_coverage(script_dir + rb_model_diam_dir + rb_cur_tree_file, trees[index][tree_counter])
        print(f"rb, diam={diameter}, tree={tree_counter}, average_posterior_rb_distance={rb_coverage}")
        print(f"TIP, diam={diameter}, tree={tree_counter}, average_posterior_rb_distance={TIP_coverage}")
