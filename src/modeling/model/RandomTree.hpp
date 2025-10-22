#include "modeling/parameters/RateMatrix.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "core/EigenSystem.hpp"
#include "core/RateEigens.hpp"
#include <filesystem>


class ConditionalLikelihood;
class TransitionProbability;

/* This classs exists to generate a validation tree where the true topology is known. When this class is initialized
   The tree topology, rateMatrix and stationary are fixed. When genNewData is called, a fasta file with the tree tip
   sequences will be outputed to XX FOLDER. Branch lengths can be changed
*/
class RandomTree {
    public:
        RandomTree(int numTaxa, int characterSequenceCount, double lambda, std::filesystem::path outputPath); // initialize a random tree and the sequence outputs.
        ~RandomTree(); // randomTree destructor
        void genNewData(); // generate new sequence data based on randomTree parameters and store them in allNodeSequences
        void changeBranchLengths(double newRateParam); // change the branch lengths by drawing new ones from an exponential with parameter newRateParam
        std::vector<Node*> getPreorderTrav(); // finds and returns the pre-order traversal of the tree as a vector of nodes
        TreeObject* getTree() {return tree;} // return the tree object
        std::vector<std::vector<int>> getAllNodeSequences(){return allNodeSequences;} // return all sequences for all nodes as a vector of vector
        int printSequences(std::string fileName); // take the sequences stored in allNodeSequences and print them outputPath/fileName
        int printTips(std::string fileName); // take the tip/leaf sequences stored in allNodeSequences and print them outputPath/fileName
        void rescaleTree(double newTreeDiameter); // rescale the treediameter to newTreeDiameter
        std::string getOutputPath() {return outputPath;} // return the output path
        void setOutputPath(std::string newOutputPath) {outputPath = newOutputPath;} // set a new output path
    private:
        void preorderDescend(std::vector<Node*>& PreOrderTrav, Node* currentNode); // helper method for getPreorderTrav
        RateMatrix* rateMatrix; // rateMatrix used to generate the Q matrix for transition probabilities
        TransitionProbability* transProb; // transition probability object used to generate transition probability matrices along branches, used to generate sequences
        TreeObject* tree; // tree object representing the topology
        std::vector<double> branchLengths; // vector of branch lengths for each branch in the tree
        double expRateParam; // governs the exponential distribution used to generate branch lengths
        int numTaxa; // number of taxa in the tree
        int characterStateCount; // how many nucleotides in each sequence
        RandomVariable& rng; // store the rng instance
        EigenSystem *eigens = new EigenSystem(4); // eigensystem used to diagonalize Q matrices and calculate transition probabilities
        RateEigen *rateEigen = new RateEigen(4); // stores the rate eigen decomposition for real eigenvalues
        ComplexRateEigen *complexRateEigen = new ComplexRateEigen(4); // stores the rate eigen decomposition for complex eigenvalues
        std::vector<std::vector<int>> allNodeSequences; // stores the sequences for all nodes in the tree, indexed by node index
        std::string outputPath; // path to output the generated sequences
};