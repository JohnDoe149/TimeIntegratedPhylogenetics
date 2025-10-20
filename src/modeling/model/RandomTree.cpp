#include "RandomTree.hpp"
#include "TransitionProbability.hpp"
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <numbers>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

/*
Start by initializing a random topology, random stationary and random rate matrix. numTaxa of 50 with a 
sequenceLength of 200 is a good start.
*/
RandomTree::RandomTree(int numTaxa, int characterStateCount, double lambda, std::filesystem::path newOutputPath): numTaxa(numTaxa), tree(new TreeObject(numTaxa)), transProb(new TransitionProbability(numTaxa)), 
      rateMatrix(new RateMatrix()), expRateParam(lambda), characterStateCount(characterStateCount), rng(RandomVariable::randomVariableInstance()),
      allNodeSequences(numTaxa * 2 - 1)
      {
    
    // Each branch length for the topology is drawn from an exponential prior
    // each node in the tree (except the root) has an ancestral branch, so use the index of the node 
    // to index into branchLengths to get the branchLength
    for(int i = 0; i < (numTaxa*2-1); i++){
        branchLengths.push_back(Probability::Gamma::rv(&rng, 1, expRateParam));
    }
    tree->setNodeNameIndex();
    tree->setBranchLength(branchLengths);
    outputPath = newOutputPath;
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
    Matrix<std::complex<double>> cDiagonalMatrix(4, 4, 0.0);
    Matrix<double> diagonalMatrix(4, 4, 0.0);
    Matrix<double> preDiag = rateMatrix->Q();
    bool isComplex = eigens->update(preDiag, *rateEigen, *complexRateEigen);
    for(int i = 0; i < 4; i++){
        cDiagonalMatrix(i, i) = complexRateEigen->ceigenvalue[i];
        diagonalMatrix(i, i) = rateEigen->eigenvalue[i];
    }

    // skip the root node, so start nodeTrav at 1
    for(int nodeTrav = 1; nodeTrav < treeTraversal.size(); nodeTrav++){
        Node* currNode = treeTraversal[nodeTrav];

        // now calculate the transition probabilties, e^(rateMatrix * branchlength) using TransitionProbability
        Matrix<double> transProbMatrix(4, 4, 0);
        if(isComplex){
            
            // we will use tempComplex as an intermediate matrix for calculations
            Matrix<std::complex<double>> tempComplex(4, 4, 0.0);

            // exponentiate tempComplex
            for(int i = 0; i < 4; i++){
                // multiply the diagonal matrix by the branch lengths before exponentiating it
                tempComplex(i,i) = std::exp(cDiagonalMatrix(i, i) * branchLengths[currNode->getIndex()]);
            }

            // get the original matrix back
            Matrix<std::complex<double>> newMatrix = (( (*complexRateEigen->cDiagLeftMatrix) * tempComplex)*(*complexRateEigen->cDiagRightMatrix));  
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    transProbMatrix(i, j) = newMatrix(i, j).real();
                }
            }
        } else {
            // exponentiate diagonalMatrix
            Matrix<double> tempMatrix(4, 4, 0.0);

            for(int i = 0; i < 4; i++){
                double temp = diagonalMatrix(i, i) * branchLengths[currNode->getIndex()];
                tempMatrix(i, i) = std::exp(temp);
            }

            // get the original matrix back
            Matrix<double> newMatrix = (( (*rateEigen->diagLeftMatrix) * tempMatrix)*(*rateEigen->diagRightMatrix));          
            transProbMatrix = newMatrix;
        }
        
        
        int ancIndex = currNode->getAncestor()->getIndex();
        std::vector<int> currNodeAncestorSeq = allNodeSequences[ancIndex]; // get the ancestor sequence

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

    // next go through the current node's neighbors and recursively call the method
    // on the current node's descendants
    std::set<Node*>& nodeNeighbors = currentNode->getNeighbors();
    for(Node* n : nodeNeighbors){
        if(n!=currentNode->getAncestor()){
            preorderDescend(preOrderTrav, n);
        }
    }
}

// a method to output the generated sequence to a file. Returns 1 if there is an error creating the file
// Outputs a file to the testing folder and if a file with the filename already exists, overwrites it.
// The parameter fileName should NOT have a format (e.g., .fasta/.txt). NOTE: this prints the sequences for all the trees
// not just the tips
int RandomTree::printSequences(std::string fileName){
    std::ofstream outFile((outputPath + "/" + fileName));
    for(int taxaIter = 0; taxaIter < allNodeSequences.size(); taxaIter++){
        std::vector<int> currentSequence = allNodeSequences[taxaIter];
        outFile << ">Taxa_number_" << taxaIter << "\n";
        for(int seqIter = 0; seqIter < currentSequence.size(); seqIter++){
            if(currentSequence[seqIter] == 0){
                outFile << "A";
            } else if (currentSequence[seqIter] == 1){
                outFile << "C";
            } else if (currentSequence[seqIter] == 2){
                outFile << "G";
            } else if (currentSequence[seqIter] == 3){
                outFile << "T";
            }   
        }
        outFile << "\n";
    }
    return 0;
}

// a method to output the generated sequences. Returns a 1 if there is an error creating the file. 
// Outputs a file to the testing folderand if a file with the filename already exists, overwrites it.
// the fileName should not have a format. NOTE: this should be used for simulation and testing as it prints 
// just the tips.
int RandomTree::printTips(std::string fileName){
    std::ofstream outFile((outputPath + "/" + fileName));
    // iterate through each tip node of the generated tree.
    std::vector<Node*> allTips = tree->getTips();
    for(int taxaIter = 0; taxaIter < allTips.size(); taxaIter++){

        // fortunately, all of the nodes have an index that corresponds to their sequence in "allNodeSequences".
        // therefore, go through each of the generated tree's tips and get their index to access their sequence.
        // this has the effect of only printing sequences that come from tip nodes.
        std::vector<int> currentSequence = allNodeSequences[allTips[taxaIter]->getIndex()];
        outFile << ">Taxa_number_" << taxaIter << "\n";
        for(int seqIter = 0; seqIter < currentSequence.size(); seqIter++){
            if(currentSequence[seqIter] == 0){
                outFile << "A";
            } else if (currentSequence[seqIter] == 1){
                outFile << "C";
            } else if (currentSequence[seqIter] == 2){
                outFile << "G";
            } else if (currentSequence[seqIter] == 3){
                outFile << "T";
            }   
        }
        outFile << "\n";
    }
    return 0;
}

// a method to rescale randomTree's randomly generated tree's diameter to the parameter newTreeDiameter. 
// The old tree diameter is found by taking the longest branch from root to tip. This is down with a DFS (preorder traversal)
// Call this method before generating new data.
void RandomTree::rescaleTree(double newTreeDiameter){
    Node* rootNode = tree->getRoot();
    std::vector<Node*> visitNode;
    visitNode.push_back(rootNode);
    std::vector<double> cumulativeBranchLengths(branchLengths);

    while(visitNode.size() != 0){

        // get the current node and add their children to the queue
        Node* currentNode = visitNode.back();
        visitNode.pop_back(); 
        std::set<Node*> neighborSet = currentNode->getNeighbors();
        neighborSet.erase(currentNode->getAncestor());
        for (Node* childNode: neighborSet){
            visitNode.push_back(childNode);
        }

        // now calculate the cumulative branch length by calling on the ancestor of the current node
        // safety check if the currentNode is root
        if(currentNode != rootNode){
            cumulativeBranchLengths[currentNode->getIndex()] += cumulativeBranchLengths[currentNode->getAncestor()->getIndex()];
        } 
    }

    // now all of the branchlengths have been counted, find the longest diameter by iterating through all of the tips
    // also keep track which tip node it is so we can find how many branches up the diameter
    std::vector<Node*> tips = tree->getTips();
    double maxDiameter = 0;
    Node* longestNode;
    for(Node* node : tips){
        double tentDiameter = cumulativeBranchLengths[node->getIndex()];
        if(tentDiameter > maxDiameter){
            maxDiameter = tentDiameter;
            longestNode = node;
        }
    }

    // calculate what we need to multiply each branch length by to achieve our desired treediameter
    double scale = newTreeDiameter/(maxDiameter);
    for(int index = 0; index < branchLengths.size(); index++){
        branchLengths[index] *= scale;
    }
    tree->setBranchLength(branchLengths);
    return;
}


RandomTree::~RandomTree(){
    delete tree;
    delete transProb;
    delete rateMatrix;
}