#include "RandomTree.hpp"
#include "TransitionProbability.hpp"
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <numbers>
#include <cmath>


/*
Start by initializing a random topology, random stationary and random rate matrix. numTaxa of 50 with a 
sequenceLength of 200 is a good start.
*/
RandomTree::RandomTree(int numTaxa, int characterStateCount) 
    : numTaxa(numTaxa), tree(new TreeObject(numTaxa)), transProb(new TransitionProbability(numTaxa)), 
      rateMatrix(new RateMatrix()), expRateParam(1.1), characterStateCount(characterStateCount), rng(RandomVariable::randomVariableInstance(100)),
      allNodeSequences(numTaxa * 2 - 1) 
      {
    
    // Each branch length for the topology is drawn from an exponential prior
    // each node in the tree (except the root) has an ancestral branch, so use the index of the node 
    // to index into branchLengths to get the branchLength
    for(int i = 0; i < (numTaxa*2-1); i++){
        branchLengths.push_back(Probability::Gamma::rv(&rng, 1, expRateParam));
    }
}

/* This method generates a new set of sequence data
*/
void RandomTree::genNewData(){

    // first get the root of the topology and generate the ancestral root sequence, we will use ATCG where
    // A = 0, C = 1, G = 1, T = 1
    Node* root = tree->getRoot();
    std::vector<double> stationary = rateMatrix->getStationary();

    // draw a number between 0 and 1 to determine the character state at each site
    std::vector<int> rootSequence(characterStateCount);
    for(int i = 0; i < characterStateCount; i++){
        double draw = rng.uniformRv();
        double cumulative = 0;
        
        // loop through the possible character states to figure out where the draw landed us
        // this is possible because stationary sums to 1 and rng.uniformRv draws between 0 and 1
        for(int state = 0; state < stationary.size(); state++){
            cumulative += stationary[state];
            if(draw < cumulative){
                rootSequence[i] = state;
                break;
            }
        }
    }
    allNodeSequences[root->getIndex()] = rootSequence;

    // start the pre-order traversal to start generating some the tip data.
    // with the ancestor sequence, calculate the transition probability matrix to get 
    std::vector<Node*> treeTraversal = getPreorderTrav();

    // skip the root node, so start nodeTrav at 1
    for(int nodeTrav = 1; nodeTrav < treeTraversal.size(); nodeTrav++){
        Node* currNode = treeTraversal[nodeTrav];

        // now calculate the transition probabilties, e^(rateMatrix * branchlength) using TransitionProbability
        Matrix<double> transProbMatrix(4, 4, 0);
        Matrix<double> preDiag = rateMatrix->Q() * branchLengths[currNode->getIndex()];
        preDiag.print();
        bool isComplex = eigens->update(preDiag, *rateEigen, *complexRateEigen);
        if(isComplex){
            Matrix<std::complex<double>> diagonalMatrix(4, 4, 0.0);
            for(int i = 0; i < 4; i++){
                diagonalMatrix(i, i) = complexRateEigen->ceigenvalue[i];
            }

            // exponentiate diagonalMatrix
            for(int i = 0; i < 4; i++){
                diagonalMatrix(i, i) = std::pow(M_E, diagonalMatrix(i, i));
            }

            // get the original matrix back
            Matrix<std::complex<double>> newMatrix = (( (*complexRateEigen->cDiagLeftMatrix) * diagonalMatrix)*(*complexRateEigen->cDiagRightMatrix));  
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    transProbMatrix(i, j) = diagonalMatrix(i, j).real();
                }
            }
        } else {
            Matrix<double> diagonalMatrix(4, 4, 0.0);
            for(int i = 0; i < 4; i++){
                diagonalMatrix(i, i) = rateEigen->eigenvalue[i];
            }
            diagonalMatrix.print();

            // exponentiate diagonalMatrix
            for(int i = 0; i < 4; i++){
                diagonalMatrix(i, i) = std::pow(M_E, diagonalMatrix(i, i));
            }
            diagonalMatrix.print();

            // get the original matrix back
            Matrix<double> newMatrix = (( (*rateEigen->diagLeftMatrix) * diagonalMatrix)*(*rateEigen->diagRightMatrix));          
            transProbMatrix = newMatrix;
        }
        
        int ancIndex = currNode->getAncestor()->getIndex();
        std::vector<int> currNodeAncestorSeq = allNodeSequences[ancIndex]; // get the ancestor sequence
        transProbMatrix.print();

        // now calculate the currentnode's sequence from the transProbMatrix
        std::vector<int> currentSequence;
        std::vector<double> transitionProbVector;
        for(int i = 0; i < characterStateCount; i++){
            transitionProbVector.clear();

            // now get the correct row from the transition probability matrix
            for(int j = 0; j < 4; j++){
                transitionProbVector.push_back(transProbMatrix(currNodeAncestorSeq[i], j));
            }
            double draw = rng.uniformRv();
            double cumulative = 0;
            
            // draw the current characterStateCount from the ancestor
            for(int state = 0; state < transitionProbVector.size(); state++){
                cumulative += transitionProbVector[state];
                if(draw < cumulative){
                    currentSequence.push_back(state);
                    break;
                }
            }
        }

        // we now have the new currentSequence, so add it to the array of sequences
        allNodeSequences[currNode->getIndex()] = currentSequence;
    }
}

// This method lets you change the expRateParam and generate new branch lengths with it.
void RandomTree::changeBranchLengths(double newRateParam){
    expRateParam = newRateParam;
    branchLengths.clear();

    for(int i = 0; i < (numTaxa*2-3); i++){
        branchLengths.push_back(Probability::Gamma::rv(&rng, 1, newRateParam));
    }
}

// this method starts a pre-order traversal
std::vector<Node*> RandomTree::getPreorderTrav(){
    std::vector<Node*> preOrderTrav;

    // call preorderDescend on the root
    preorderDescend(preOrderTrav, tree->getRoot());
    return preOrderTrav;
}

// the recursive part of the method to get the pre-order traversal
void RandomTree::preorderDescend(std::vector<Node*>& preOrderTrav, Node* currentNode){
    if(currentNode == nullptr){
        return;
    }

    // pushback the current node
    preOrderTrav.push_back(currentNode);
    std::set<Node*>& nodeNeighbors = currentNode->getNeighbors();
    for(Node* n : nodeNeighbors){
        if(n!=currentNode->getAncestor()){
            preorderDescend(preOrderTrav, n);
        }
    }
}

RandomTree::~RandomTree(){
    delete tree;
    delete transProb;
    delete rateMatrix;
}