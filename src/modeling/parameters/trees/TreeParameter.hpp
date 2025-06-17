#ifndef TREE_PARAMETER_HPP
#define TREE_PARAMETER_HPP
#include "modeling/parameters/Parameter.hpp"
#include "TreeObject.hpp"
#include <string>

class TreeParameter : public Parameter{
    public:
        TreeParameter(void)=delete;
        TreeParameter(Alignment* aln, std::string newick, double lambda);
        ~TreeParameter();
        TreeObject* getTree(){return trees[0];}

        // For making topology changes with NNI
        void acceptTreeMove();
        double updateTreeMove();

        // for making length changes t
        void acceptTreeGamma();
        double updateTreeGamma();

        void reject();
        void tune();
        double lnPrior();

        std::string writeNewick() {return trees[0]->getNewick();}

        int branchCount;
        int branchAcceptCount;
        int treeCount; 
        int treeAcceptCount;
    private:
        bool fixedTree;
        int moveChoice;
        double branchDelta;
        double treeAlpha;
        double lambda;
        double currentPrior;
        double oldPrior;
        TreeObject* trees[2];
};

#endif