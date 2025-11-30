#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <string>

class TreeParameter : public Parameter{
    public:
        TreeParameter(void)=delete;
        TreeParameter(Alignment* aln); // initializes a treeObject based on the alignment and branch gammas
        ~TreeParameter(); // delete the stored tree objects
        TreeObject* getTree(){return trees[0];} // return the current/active tree object
        double updateTreeMove(); // Performs one NNI move on the current tree
        double updateTreeGamma(); // Rescales the shape or rate parameter of a randomly selected branch
        void accept(); // accept the current tree and its parameters
        void reject(); // reject the current tree and its parameters and revert to the previous accepted tree
        void tune(); // tune the proposal step-sizes for branch rate and shape parameters
        double lnPrior(){return currentPrior;} // return the prior for the branch gamma parameters
        std::string writeNewick() {return trees[0]->getNewick();} // return the newick string of the current tree
    private:
        // store the current and last accepted prior
        double currentPrior;
        double oldPrior;

        // store two tree objects, one for current proposals and one for storing the last accepted tree
        TreeObject* trees[2];

        // branch gamma prior parameters
        double shapePriorRate;
        double ratePriorRate;
        double shapePriorShape;
        double ratePriorShape;

        // proposal step-sizes
        double shapeDelta;
        double rateDelta;

        // tracking acceptance counts for burn-in tuning
        int moveChoice;
        int shapeCount;
        int shapeAcceptCount;
        int rateCount;
        int rateAcceptCount;
};

#endif