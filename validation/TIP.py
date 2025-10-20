import subprocess
import os

diams = ["0.250000", "0.500000", "1.000000", "2.000000", "3.000000", "4.000000"]
samples_per_partition = 5
script_dir = os.path.dirname(os.path.abspath(__file__))
build_dir = os.path.join(script_dir,"..","build")
file_dir = os.path.join(script_dir,"..", "validation")
build_dir = os.path.abspath(build_dir)
file_dir = os.path.abspath(file_dir)
print(build_dir)

# loop through the different partitions of data that are partitioned by tree length
for diam_index, diam in enumerate(diams):
    for sample_index in range(samples_per_partition):
        input_file = f"/TIP_testing_data/taxa50diam{diam}siteCount800/{sample_index}tipSequence.fasta"
        output_trace_file = f"/workspaces/JohnPhylo/validation/output/taxa50diam{diam}siteCount800/TIP{sample_index}tipSequence.log"
        output_tree_file = f"/workspaces/JohnPhylo/validation/output/taxa50diam{diam}siteCount800/TIP{sample_index}tipSequence.tree"

        # now run the subprocess
        subprocess.run(["./main", "-nexus", file_dir+input_file, "-treeOut", output_tree_file, "-mcmcOut", output_trace_file], cwd=build_dir)
