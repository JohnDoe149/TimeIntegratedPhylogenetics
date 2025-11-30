#ifndef CONDITIONAL_LIKELIHOOD_HPP
#define CONDITIONAL_LIKELIHOOD_HPP
#include <iostream>

class Alignment;
class Node;

 class ConditionalLikelihood{
    public:
        ConditionalLikelihood(void) = delete; // remove default constructor
        ConditionalLikelihood(Alignment* aln, int nN); // aln stores alignment, used to get number of characters in sequence, nN determines number of nodes in tree. Determines how much memory allocated for CL
        ~ConditionalLikelihood(); // destructor
        double* operator()(int n); // returns a pointer to the conditional likelihoods for node with n
    private:
        double* condLikelihoods; // pointer to all conditional likelihoods stored as a flat array
        int numChar; // number of characters in the alignment for each node
        int numNodes; // number of nodes in the tree
        int stateSpace; // number of states in the model (in this case 4, for nucleotides)
};

#endif