#ifndef CONDITIONAL_LIKELIHOOD_HPP
#define CONDITIONAL_LIKELIHOOD_HPP
#include <iostream>

class Alignment;
class Node;

 class ConditionalLikelihood{
    public:
        ConditionalLikelihood(void) = delete;
        ConditionalLikelihood(Alignment* aln, int nN);
        ~ConditionalLikelihood();
        double* operator()(int n);
    private:
        double* condLikelihoods;
        int numChar;
        int numNodes;
        int stateSpace;
};

#endif