# Time Integrate Phylogenetics

## Authors and Acknowledgement 
This project is built on a legacy phylogenetic codebase originally designed for complex Bayesian phylogenetic inference, odels. The first contributor to the legacy code base is John Huelsenbeck and colleagues. They wrote the code for the EigenSystem.cpp/hpp, Math.cpp/hpp and Matrix.hpp that I leverage to do eigendecomposition. The second contributor built the majority of the code required for phylogenetic inference, Wesley C Demontigny (my mentor). He stripped down an older, more specialized inference model to allow me to freely make use of and modify classes like TreeParameter and RateMatrix classes. These provided classes would serve as the skeleton for my model. I'd like to thank Wesley Demontigny for being a mentor to me and writing the majority of phylogenetic Bayesian inference code. Without the past contributions of my mentor, this would never have been possible.

## Introduction
The goal was to develop a phylogenetic Bayesian inference model that would be more robust under pathological data. An example of incorrect data would be errors in genome assembly or mislabeling of data. Wesley guided me on the design and it was decided that the model would take advantage of the fact that if a branch length has a gamma-prior, then the transition probability would follow the formula 
P(Q, alpha, beta) = (I - 1/beta * Q)^ -alpha. So In order to make the model more robust, each branch of the tree would have its length drawn from its own gamma distribution. In other words, each branch would have its own gamma parameters associated with it. The original idea was that this extra "layer" and the fact that the gamma distribution gets "integrated out" in the original derivation for the transition probability would make the model more resistant/robust to incorrect data. By having the distribution "integrated out", the model essentially considers all possible branch lengths that are possible from the distribution. 
The novelty of this model is that each branch is independent for each other and has a gamma distribution for its branch lengths. These gamma distribution's are integrated out to consider all possible branch lengths when evaluating the model and sampling from the posterior. In order to implement this novel model, a number of improvements to the original codebase had to be made. This includes stripping down extraneous features, writing new methods, implementing robust scenario testing and extensive documentation.

## My Contributions (John Du)
- modified rateEigen.hpp and EigenSystem.cpp/hpp to store and perform eigendecomposition for calculations
- Implemented a Blocked Gibbs-sampler to improve chain convergence to true topology, improved accuracy by 40% over unblocked scheme
- Implemented GTR-substitution model to replace base JC-substitution model, including initialization, priors, proposal and tuning for the transition rates.
- Implemented a new type of tree topology proposal called NNI (Nearest neighbor interchange) in TreeParameter.cpp
- Implemented branch gamma parameters to replace traditional branch lengths, which includes initialization, priors, proposals and tuning
- Implemented Huelsenbeck's transition probability formula and integrated it into Felsenstein's pruning algorithm for conditional likelihood
- Wrote extensive documentation for every modified method and class (code was previously scarcely documented)
- Reduced redundancy and improved readability in code base by trimming down
- Implemented robust testing for modified and new features such as eigen-decomposition
- Tuned hyperparameters such as move weights and gibb's iterations.
- Wrote a custom class to generate MSA (Multiple Sequence Alignments) from known trees for validation to see if model can deduce the true tree given data.
- Wrote custom scripts to run this model and revBayes GTR on custom MSAs generated from knwon true trees
- Wrote validation scripts leveraging libraries like DendroPy to evaluate model performance


## Challenges 
Part of the challenge was understanding the inference framework; while I understood the theory behind MCMC, actually understanding how it can be implemented in a practical manner enough to build something new is a different kind of task. As a result, I grew better at reading code I did not write and became more aware of the quirks in my own code. 

When I first started working on this project, I had no experience with C++, only some familiarity with C. So I was forced to adapt and quickly pick up on the nuances that seperate C from C++. C++'s object oriented patterns, memory management and differences in standard library usage were all things I became familiar with as I worked on this project.

Computation time quickly became an issue due to the nature of the method employed. Eigendecomposition is expensive computationally, so while a typical inference model also performs eigendecomposition, typical models can reuse the eigenvalues so long as the Q matrix does not change. However, due to how my model calculates transition probabilities, eigendecomposition must be performed if either the branch parameters or Q matrix changes. Branch parameters are something that must be updated frequently, so this is hurts computation speed and this combined with other factors could result in relatively small runs (100,000 iterations on a small tree) taking upwards of several hours.

## Results
Ultimately, the hope was for the model to perform better on datasets generated from longer tree diameters (less informative). So I conducted 25 trial runs that pitted my model against a standard inference model (revBayes' GTR) across 5 treeDiameters ranging from 0.25 to 4.00. My model would perform similarly to revBaye's GTR, but with considerable worse run time. The worse runtime can be attributed to factors like the nature of Huelsenbeck's transition probability formula, greatly increasing the number of eigendecompositions required. While this is done with relatively few data points, given that my model is rarely able to outcompete revBayes and the fact that it demands far more resources computationally, the gamma-parameterized branch length no-common-mechanism model is not worth using. It produces about the same performance but with a significant tradeoff in terms of run time. 

## Structure and Navigation

Inference models like this one output tree and trace files at the end of a run, the location and name of such files is determined via terminal commandline but can also be preset in "Settings.cpp"
Validation scripts to analyze tree files and automate model runs can be found in the "validation" folder. That is also where generated simulation MSAs are stored.
"src" contains all of the core files used to run the model. The folder "core" with src contains all the C++ files used primarily for math applications such as EigenSystem and Matrix. The folder "modeling" contains three folders.
"analysis" folder contains the Mcmc.cpp which handles the overall flow of mcmc. "model" folder contain the C++ files used to sample from the posterior such as ConditionalLikelihood and TransitionProbability. "parameter" folder has all of the parameters related files like the treeParameter and rateMatrix.
Also within "src" is the main.cpp file. This is where all the model is initialized. 

## Support
If you have any questions regarding the model, send an email to: johndu@terpmail.umd.edu
