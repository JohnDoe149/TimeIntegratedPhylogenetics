# Time Integrated Phylogenetic Bayesian Inference Model (TIP-BIM)

## Authors and Acknowledgement 
This project is built on a legacy phylogenetic codebase originally designed for complex Bayesian Models. The first contributor to the legacy code base is John Huelsenbeck (who is the mentor of the second contributor) and colleagues. They wrote the code for the EigenSystem.cpp/hpp, Math.cpp/hpp and Matrix.hpp that allow for an easy way to do tedious and complex matrix transformations. The second contributor built the majority of the code required for
phylogenetic inference, Wesley C Demontigny (my mentor). The Alignment, Probability, RateEigens, Mcmc, ConditionalLikelihood, Model, TransitionProbability, Node, TreeObject, TreeParameter and RateMatrix classes are all key contributions he made that served as the skeleton for my own model. I'd like to thank Prof. Huelsenbeck for the robust diagonalization methods offered by EigenSystem. I'd like to thank Wesley Demontigny for being a mentor to me and writing the majority of phylogenetic Bayesian inference code. Without the past contributions of my mentor, this would never have been possible.

## Introduction
The goal was to develop a phylogenetic Bayesian inference model that would be more robust under pathological/incorrect data. An example of incorrect data would be errors in genome assembly or mislabeling of data. Wesley Demontigny guided me on the design and it was decided that the model would take advantage of the fact that if a branch length has a gamma-prior, then the transition probability would follow the formula 
P(Q, alpha, beta) = (I - 1/beta * Q)^ -alpha. So In order to make the model more robust, each branch of the tree would have its length drawn from its own gamma distribution. In other words, each branch would have its own gamma parameters associated with it. The original idea was that this extra "layer" and the fact that the gamma distribution gets "integrated out" in the original derivation for the transition probability would make the model more resistant/robust to incorrect data. By having the distribution "integrated out", the model essentially considers all possible branch lengths that are possible from the distribution. 
The novelty of TIP model is that each branch is independent for each other and has a gamma distribution for its branch lengths. These gamma distribution's are integrated out to consider all possible branch lengths when evaluating the model and sampling from the posterior. In order to implement this novel model, a number of improvements to the original codebase had to be made. This includes stripping down extraneous features, writing new methods, implementing robust scenario testing and extensive documentation.

## My Contributions (John Du)
- modified rateEigen.hpp and EigenSystem.cpp/hpp to easily diagonalize a matrix with complex eigenvalues
- Implemented a new type of tree topology proposal called NNI (Nearest neighbor interchange) in TreeParameter.cpp
- Implemented a alpha and beta parameter for each branch's gamma distribution
- Wrote updating, tuning and a Gibbs sampler for each branch's alpha and beta parameters
- Implemented the integrated form of the Transition Probability matrix calculation for branch's with gamma priors
- Wrote extensive documentation for every modified function and file
- Implemented robust testing for modified and new features


## Challenges 
While Huelsenbeck's contributions were well documented, I could unfortunately not say the same thing about Demontigny's contributions. One of the challenges I faced was interpreting undocumented code and understanding the design quirks. The accept-reject memory buffer design was a design I was unfamiliar with and due to the lack of documentation I was unable to understand for some time. This resulted in an extensive effort to reverse engineer the purpose and functionality various aspects of the existing codebase. The location of where key variables are located, seemingly random static class variables that were actually vital and foreign design patterns were all challenges I ran into.
When I initially began working on this model, I had no experience working in C++. I only knew C at the time 
and was forced to adapt and quickly pick up on the nuances that seperate C from C++. C++'s object oriented patterns, memory management differences and differences in standard library usage all served as hurdles that I had to overcome.

## Structure and Navigation
From the root directory, there are two folders of interest. 
"res" contains the "analysis.log" file which can be read by applications such as Tracer to view the trace plots of the model runs. "trees.trees" contains the tree file that can be read by applications such FigTree to view the tree that the model produced.
"src" contains all of the core files used to run the model. The folder "core" with src contains all the C++ files used primarily for math applications such as EigenSystem and Matrix. The folder "modeling" contains three folders.
"analysis" folder contains the Mcmc.cpp which handles the overall flow of mcmc. "model" folder contain the C++ files used to sample from the posterior such as ConditionalLikelihood and TransitionProbability. "parameter" folder has all of the parameters related files like the treeParameter and rateMatrix.
Also within "src" is the main.cpp file. This is the file that is run to run the model. 

## Usage

## Support
If you have any questions regarding the model, send an email to: johndu@terpmail.umd.edu
