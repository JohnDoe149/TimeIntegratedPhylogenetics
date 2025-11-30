#ifndef PHYLO_CTMC_HPP
#define PHYLO_CTMC_HPP
#include "modeling/parameters/trees/TreeParameter.hpp"
#include <taskflow/taskflow.hpp>
#include "core/Alignment.hpp"

class ConditionalLikelihood;
class TransitionProbability;
class RandomVariable;
class RateMatrix;
class Settings;

class Model {
    public:
        Model(void) = delete; // remove default constructor
        Model(Settings s, Alignment* a, TreeParameter* t, RateMatrix* m); // initialize model with alignment, treeParameter and rateMatrix as parameters
        ~Model(); // destructor

        double lnLikelihood() {return currentLikelihood;} // return lnLikelihood for posterior calculations
        double lnPrior(); // return lnPrior as the sum of lnPrior from its parameters like Tree and RateMatrix

        void regenerateLikelihood(); // regenerate the likelihood to update currentLikelihood

        int getNumTaxa(){return aln->getNumTaxa();} // return numTaxa
        int getNumChar(){return numChar;} // return numChar
        int getNumNodes(){return numNodes;} // return numNodes

        TransitionProbability* getTransitionProbability() { return transProb; } 
        ConditionalLikelihood* getConditionalLikelihood() { return postOrder; }

        void accept();
        void reject();
        void tuneMoves();

        std::string tabularHeader();
        std::string tabularOut(int i);
        std::string treeHeader();
        std::string treeOut(int i);
    protected:
        double oldLikelihood;
        double currentLikelihood;
    private:
        int stateSpace;
        int numChar;
        int numNodes;
        bool* activeTP;
        tf::Executor executor;
        RateMatrix* rateMatrix;
        Alignment* aln;
        ConditionalLikelihood* postOrder;
        TransitionProbability* transProb;
        double* rescaling;
        TreeParameter* tree;
};

#endif