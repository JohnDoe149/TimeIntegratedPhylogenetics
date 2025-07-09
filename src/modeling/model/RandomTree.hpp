#include "modeling/parameters/RateMatrix.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/EigenSystem.hpp"
#include "core/RateEigens.hpp"


class ConditionalLikelihood;
class TransitionProbability;

/* This classs exists to generate a validation tree where the true topology is known. When this class is initialized
   The tree topology, rateMatrix and stationary are fixed. When genNewData is called, a fasta file with the tree tip
   sequences will be outputed to XX FOLDER. Branch lengths can be changed
*/
class RandomTree {
    public:
        RandomTree(int numTaxa, int characterSequenceCount);
        ~RandomTree();
        void genNewData();
        void changeBranchLengths(double newRateParam);
        std::vector<Node*> getPreorderTrav();
        TreeObject* getTree() {return tree;}
        std::vector<std::vector<int>> getAllNodeSequences(){return allNodeSequences;}
    private:
        void preorderDescend(std::vector<Node*>& PreOrderTrav, Node* currentNode);
        RateMatrix* rateMatrix;
        TransitionProbability* transProb;
        TreeObject* tree;
        std::vector<double> branchLengths;
        double expRateParam;
        int numTaxa;
        int characterStateCount;
        RandomVariable& rng;
        EigenSystem *eigens = new EigenSystem(4);
        RateEigen *rateEigen = new RateEigen(4);
        ComplexRateEigen *complexRateEigen = new ComplexRateEigen(4);
        std::vector<std::vector<int>> allNodeSequences;
};