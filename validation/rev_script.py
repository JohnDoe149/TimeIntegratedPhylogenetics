import subprocess
import tempfile
import os

def rev_run(diam, sample_index):
    input_file = f"TIP_testing_data/taxa50diam{diam}siteCount800/{sample_index}tipSequence.fasta"
    output_trace_file = f"rb_output/taxa50diam{diam}siteCount800/rb{sample_index}tipSequence.log"
    output_tree_file = f"rb_output/taxa50diam{diam}siteCount800/rb{sample_index}tipSequence.tree"
    rev_script_template = f"""
    # Load data
    data <- readDiscreteCharacterData(file="{input_file}")
    moves = VectorMoves()

    num_sites <- data.nchar()
    data_taxa <- data.taxa() 
    num_taxa <- data.ntaxa()
    n_branches <- 2 * num_taxa - 3
    tree_topology ~ dnUniformTopologyBranchLength(data_taxa, dnExponential(10.0))

    # Exchange rates
    er ~ dnDirichlet(v(1,1,1,1,1,1))

    # Base frequencies
    pi ~ dnDirichlet(v(1,1,1,1))

    # Create a GTR rate matrix
    Q := fnGTR(er, pi)

    phySeq ~ dnPhyloCTMC(tree=tree_topology, Q=Q, type="DNA")
    phySeq.clamp(data)

    # Set up proposals, stationary, rates, branch lengths, and topology
    moves.append(mvSimplexElementScale(er, alpha=10.0, tune=TRUE, weight=3.0))
    moves.append(mvSimplexElementScale(pi, alpha=10.0, tune=TRUE, weight=3.0))
    moves.append(mvNNI(tree_topology, weight=8.0))
    moves.append( mvBranchLengthScale(tree_topology, weight=8.0, tune=true) )

    # Monitors
    monitors = VectorMonitors()
    monitors.append(mnModel(filename="{output_trace_file}", printgen=10, exclude=["F"]))
    monitors.append(mnFile(filename="{output_tree_file}", printgen=10, tree_topology))

    # Run MCMC
    mymodel = model(tree_topology)
    mymcmc = mcmc(mymodel, moves, monitors)
    mymcmc.burnin(generations=10000, tuningInterval=100)
    mymcmc.run(generations=50000)
    q()
    """

    # now create a temporary runPy.Rev file with the script template
    with tempfile.NamedTemporaryFile(mode= 'w+t', delete=False) as tmp:
        print(f"running on {input_file}")
        tmp.write(rev_script_template)
        tmp.flush()
        # now run the subprocess
        subprocess.run(["./rb", tmp.name])

    os.remove(tmp.name)

diams = ["0.250000", "0.500000", "1.000000", "2.000000", "3.000000", "4.000000"]
samples_per_partition = 5

# loop through the different partitions of data that are partitioned by tree length
for diam_index, diam in enumerate(diams):
    for sample_index in range(samples_per_partition):
        rev_run(diam, sample_index)

