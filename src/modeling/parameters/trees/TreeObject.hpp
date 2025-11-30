#ifndef TREE_OBJECT_HPP
#define TREE_OBJECT_HPP
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

class Node;
class Alignment;
class RandomVariable;

class TreeObject {

    public:
                            TreeObject(void) = delete; // no default constructor
                            TreeObject(int nt); // Generate a random tree with nt taxa
                            TreeObject(Alignment* aln); // generate a random tree connected to an alignment aln
                            TreeObject(const TreeObject& t); // Copy constructor
                           ~TreeObject(void); // destructor
        TreeObject&         operator=(const TreeObject& rhs); // Deep copy assignment operator

        std::vector<double> getGammaParams(Node* n) const; // return the vector of gamma parameters for a given node
        std::map<Node*, std::vector<double>> getGammaMap(); // return the map of node to gamma parameters
        std::vector<std::vector<double>> getGammas(); // return all the gamma parameters as a vector of vectors
        void                setGammaDist(Node* n, double alpha, double beta); // set the gamma distribution parameters for a given node
        
        Node*               getNodeWithIndex(int index); // returns the node with the given index
        std::string         getNewick() const; // return the newick string of the tree
        int                 getNumTaxa(){return numTaxa;} // return the number of taxa
        int                 getNumNodes(){return nodes.size();}  // return the number of nodes
        std::vector<Node*>& getPostOrderSeq() {return postOrderSeq;} // return the post-order sequence of nodes as a vector
        Node*               getRoot() {return root;} // return the root node
        std::vector<Node*>  getTips(); // returns all the tip/leaf nodes as a vector
        void                initPostOrder(void); // reinitialize the post-order sequence, useful in case of topology change
        void                passDown(Node* p, std::vector<Node*>& vec); //recursive helper method for initPostOrder
        void                updateAll(); // sets all nodes as needing transition probability updates
        void                setNodeNameIndex(); // sets all node names to their index for debugging
        void                setBranchLength(std::vector<double> newbranchLengths){branchLengths = newbranchLengths;} // takes input vector of branch lengths and replaces current branch lengths
    private:
        std::vector<double> branchLengths; // stores all the branch lengths for the tree
        std::vector<Node*>  nodes; // vector of all nodes in the tree
        int                 numTaxa; // number of taxa in the tree
        std::vector<Node*>  postOrderSeq; // vector of nodes in post-order traversal, call initPostOrder if topology changes
        Node*               root; // pointer to the root node of the tree

        Node*               addNode(void); // helper method used during tree construction
        std::map<Node*, std::vector<double>> branchGamma; // map from node to its branch gamma parameters (shape, rate)
        void                clone(const TreeObject& t); // helper method for deep copy
        void                deleteAllNodes(); // helper method to delete all nodes
        void                writeNode(Node* p, std::stringstream& strm) const; // recursive helper method for getNewick
};

#endif
