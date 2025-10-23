#include "ConditionalLikelihood.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"

// stores the conditional likelihoods for all nodes in the tree, called by model to initialize and update Cls
ConditionalLikelihood::ConditionalLikelihood(Alignment* aln, int nN) : numNodes(nN), stateSpace(4) {
    numChar = aln->getNumChar();
    int width = numNodes*numChar*stateSpace; 
    condLikelihoods = new double[width];

    for(int i = 0; i < width; i++){
        condLikelihoods[i] = 0.0;
    }

    for(int index = 0; index < aln->getNumTaxa(); index++){
        double* p = (*this)(index);
        for(int i = 0; i < numChar; i++){
            unsigned long long int state = aln->getMatrix()[index][i];
            unsigned long long int mask = 1;
            bool assigned = false;
            for(int j = 0; j < stateSpace; j++) {
                if((mask & state) != 0){
                    *p = 1.0;
                    assigned = true;
                }
                mask <<= 1;
                p++;
            }

            if(assigned == false){
                Msg::error("Never assigned a conditional value at (" + std::to_string(index) + ", " + std::to_string(i) + ")! This has state value " + std::to_string(state));
            }
        }
    }
}

// destructor
ConditionalLikelihood::~ConditionalLikelihood(){
    delete [] condLikelihoods;
}

// returns a pointer to the conditional likelihoods for node n
double* ConditionalLikelihood::operator()(int n){
    return condLikelihoods + n*numChar*stateSpace;
}
