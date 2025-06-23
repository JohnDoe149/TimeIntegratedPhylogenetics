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
        Model(void) = delete;
        Model(Settings s, Alignment* a, TreeParameter* t, RateMatrix* m);
        ~Model();

        double lnLikelihood() {return currentLikelihood;}
        double lnPrior();

        void regenerateLikelihood();

        int getNumTaxa(){return aln->getNumTaxa();}
        int getNumChar(){return numChar;}
        int getNumNodes(){return numNodes;}

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
        bool* activeCL;
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