import subprocess
import tempfile
import os

diams = ["0.250000", "0.500000", "1.000000", "2.000000", "3.000000", "4.000000"]
samples_per_partition = 5
script_dir = os.path.dirname(os.path.abspath(__file__))


# RUN THIS IN SAME DIRECTORY AS THE seq-gen executable
for diam_index, diam in enumerate(diams):
    for sample_index in range(samples_per_partition):
        input_file = f"TIP_testing_data/taxa50diam{diam}siteCount800/{sample_index}newick.tree"
        input_file = os.path.join(script_dir, input_file)
        output_sample_data_file = f"TIP_testing_data/taxa50diam{diam}siteCount800/{sample_index}tipSequence.fasta"
        output_sample_data_file = os.path.join(script_dir, output_sample_data_file)
        rev_script_template = f"""
        ./seq-gen -mGTR -l800 -n1 < {input_file} > {output_sample_data_file}
        """
        subprocess.run(rev_script_template, shell=True)
