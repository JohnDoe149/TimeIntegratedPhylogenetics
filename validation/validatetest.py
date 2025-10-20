import os
import re
from io import StringIO
from Bio import Phylo
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
tree_filename = f"/TIP_testing_data/taxa50diam0.250000siteCount800/0newick.tree"
tree = ""
with open(script_dir + tree_filename, "r") as file:

    # each file only contains one tree
    file_tree = file.read().strip()
    tree = file_tree


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
    if rf_distance == 0.0:
        return 1
    else: 
        return 0

# method to find coverage for a given tree file
def calc_coverage(treefile_path, true_tree_newick_string):
    with open(treefile_path, 'r') as tree_file:
        content = tree_file.read()

    newick_pattern = re.compile(r'\([^;]*;\s*')
    trees = newick_pattern.findall(content)

    count = []
    for cur_iter, tree in enumerate(trees):
        rf_dist = compare_tree(norm_newick(tree), norm_newick(true_tree_newick_string))
        if(rf_dist==2):
            print("close, curiter: " + str(cur_iter+1)) # do not forget that tree files are 1 indexed and python is 0 indexed
        count.append(rf_dist)
    return count

# helper method to convert Taxa_number_x to Taxax to simplify things. rb naming convention
# can be slightly different also remove annotation
def norm_newick(newick_string): 
    new_string = re.sub(r'Taxa_number_(\d+)', r'Taxa\1', newick_string)
    return re.sub(r'\[.*?\]', '', new_string)

# run a script to calculate coverage
script_dir = os.path.dirname(os.path.abspath(__file__))
rb_coverage = calc_coverage(script_dir + "/mixTree.tree", tree)
rb_coverage.sort()
rb_coverage.remove(rb_coverage[0])
plt.hist(rb_coverage, bins=20, edgecolor='black')
plt.show()
print(rb_coverage[0])


