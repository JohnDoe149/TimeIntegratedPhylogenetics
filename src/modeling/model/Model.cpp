#include "Model.hpp"
#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Msg.hpp"
#include "ConditionalLikelihood.hpp"
#include "TransitionProbability.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/RateMatrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <unordered_map>
#include <thread>
//#include <chrono>

Model::Model(Settings s, Alignment* a, TreeParameter* t, RateMatrix* m) : 
            aln(a), tree(t), rateMatrix(m), oldLikelihood(0.0), currentLikelihood(0.0), numChar(0) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    TreeObject* activeT = tree->getTree();
    stateSpace = 4;
    numChar = aln->getNumChar();

    if(aln->getNumTaxa() != activeT->getNumTaxa())
        Msg::error("Expected " + std::to_string(aln->getNumTaxa()) + 
        "taxa in the tree, but found only " + std::to_string(activeT->getNumTaxa()));
    
    numNodes = tree->getTree()->getNumNodes();

    int flagWidths = 2 * numNodes;
    activeTP = new bool[flagWidths];
    for(int i = 0; i < flagWidths; i++){
        activeTP[i] = false;
    }
    postOrder = new ConditionalLikelihood(aln, numNodes);
    transProb = new TransitionProbability(numNodes);

    int rescaleWidth = numNodes*numChar;
    rescaling = new double[rescaleWidth * 2];
    std::fill(rescaling, rescaling + rescaleWidth * 2, 0.0);

    activeT->updateAll();
}

// destructor
Model::~Model(){
    delete postOrder;
    delete transProb;
    delete [] rescaling;
    delete [] activeTP;
}

// accept the current model state by prompting each parameter such as tree and rateMatrix accept if they have changed
void Model::accept() {
    oldLikelihood = currentLikelihood; // store accepted oldLikelihood

    for(int i = 0; i < numNodes; i++){
        activeTP[i + numNodes] = activeTP[i];
    }

    std::memcpy(rescaling + numNodes*numChar, rescaling, numNodes*numChar);

    if(tree->isDirty()){
        tree->accept();
        tree->clean();
    }
    if(rateMatrix->isDirty()){
        rateMatrix->accept();
        rateMatrix->clean();
    }

    transProb->updateQ(rateMatrix->Q());
}

// restore the previously accepted likelihood and prompt model parameters to reject if they have changed.
void Model::reject() {
    currentLikelihood = oldLikelihood;

    for(int i = 0; i < numNodes; i++){
        activeTP[i] = activeTP[i + numNodes];
    }

    std::memcpy(rescaling, rescaling + numNodes*numChar, numNodes*numChar);

    if(tree->isDirty()){
        tree->reject();
        tree->clean();
    }
    if(rateMatrix->isDirty()){
        rateMatrix->reject();
        rateMatrix->clean();
    }

    transProb->updateQ(rateMatrix->Q());
}

double Model::lnPrior(){
    return tree->lnPrior() + rateMatrix->lnPrior();
}

// This method is called mainly when certain model parameters have changed (such as the rate matrix) to recalculate
// the ConditionalLikelihood of the tree. The calculation relies on having the transition probability of each node,
// So this method also checks each node that needs to have its transition probability recalculated.
void Model::regenerateLikelihood(){
    TreeObject* activeT = tree->getTree(); 

    const std::vector<Node*> poSeq = activeT->getPostOrderSeq();

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if(rateMatrix->isDirty()){
        activeT->updateAll(); // updateAll sets all nodes are not the tip to have their CL updated and all nodes to have TP updated
        transProb->updateQ(rateMatrix->Q());
    }

    //std::chrono::steady_clock::time_point rateTime = std::chrono::steady_clock::now();
    //std::cout << "Rate computation was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(rateTime - begin).count() << "[milliseconds]" << std::endl;

    // for each node, check if they need a Transition Probability update, if so recalculate their transition
    // probability and find the node's activeIndex to determine which buffer to use.
    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        if(n->getNeedsTPUpdate() == true){
            if(n != activeT->getRoot()) {
                std::vector<double> gammaVec = activeT->getGammaParams(n);
                activeTP[nIndex] ^= true;
                bool activeIndex = activeTP[nIndex];
                transProb->setProbs(activeIndex, 0, nIndex, gammaVec[0], gammaVec[1]);
            }
            n->setNeedsTPUpdate(false);
        }
    }

    //std::chrono::steady_clock::time_point probsTime = std::chrono::steady_clock::now();
    //std::cout << "Probs computation was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(probsTime - rateTime).count() << "[milliseconds]" << std::endl;

    // CL likelihood calculation, split the alignment into 100 site chunks to parallelize
    tf::Taskflow phyloTaskflow;
    const auto processor_count = std::thread::hardware_concurrency();
    int chunkSize = (int)std::ceil((double)numChar/processor_count);
    for(int range = 0; range < processor_count; range++){ 
        int start = range * chunkSize; 
        int end = start + chunkSize-1;
        end = std::min(end, numChar-1);

        phyloTaskflow.emplace([this, &poSeq, start, end](){ //define thread behavior
            int currentChunkSize = end - start + 1;
            for(Node* n : poSeq){ // iterate through each node in post order to go from tips to the root
                int nIndex = n->getIndex();
                if(n->getIsTip() == false){
                    double* pNN = (*postOrder)(nIndex) + start * stateSpace; // get the array that stores the CLs for each site for the specific node
                    std::fill(pNN, pNN + (currentChunkSize * stateSpace), 1.0); // fill the array with 1.0s the clear it ot

                    std::set<Node*>& nNeighbors = n->getNeighbors();
                    for(Node* d : nNeighbors){ // once we have a node that needs to be updated, pull from its children to calculate its new CL for each site
                        if(d != n->getAncestor()){
                            int dIndex = d->getIndex();
                            double* pN = pNN;
                            double* pD = (*postOrder)(dIndex) + start * stateSpace;
                            const Matrix<double>& P = (*transProb)(activeTP[dIndex], 0, dIndex);

                            // go through each site in the chunk
                            for(int c = 0; c < currentChunkSize; c++){
                                for(int i = 0; i < stateSpace; i++){ //state space is ATCG
                                    double sum = 0.0;
                                    for(int j = 0; j < stateSpace; j++){
                                        sum += P(i, j) * pD[j];
                                    }
                                    (*pN) *= sum;
                                    pN++;
                                }
                                pD+=stateSpace;
                            }
                        }
                    }
                    double* rescalePointer = rescaling + (numChar * nIndex) + start;
                    std::fill(rescalePointer, rescalePointer + currentChunkSize, 0.0);

                    for(int c = 0; c < currentChunkSize; c++){
                        double max = *pNN;
                        pNN++;
                        for(int i = 1; i < stateSpace; i++){
                            if(*pNN > max)
                                max = *pNN;
                            pNN++;
                        }
                        pNN -= stateSpace;
                        for(int i = 0; i < stateSpace; i++){
                            *pNN /= max;
                            pNN++;
                        }
                        *rescalePointer = std::log(max);
                        rescalePointer++;
                    }
                }
            }
        });
    }
    executor.run(phyloTaskflow).wait();


    int rIndex = activeT->getRoot()->getIndex();
    double* pR = (*postOrder)(rIndex);
    std::vector<double> f = rateMatrix->getStationary();
    double lnL = 0.0;

    for(int c = 0; c < numChar; c++){
        double like = 0.0;
        for(int i = 0; i < stateSpace; i++){
            like += pR[i]*f[i];
        }

        lnL += std::log(like);
        pR += stateSpace;
    }

    double* rescalePointer = rescaling;
    for(int i = 0, len = numNodes * numChar; i < len; i++){
        lnL += *rescalePointer;
        rescalePointer++;
    }

    currentLikelihood = lnL;

    //std::chrono::steady_clock::time_point pruneTime = std::chrono::steady_clock::now();
    //std::cout << "Pruning was completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(pruneTime - probsTime).count() << "[milliseconds]" << std::endl;
}

// This method called this current object's tree and rateMatrix and asks them to tune.
void Model::tuneMoves(){
    tree->tune();
    rateMatrix->tune();
}

// This method when called returns a string to serve as a  header that contains all of the parameters used
// in the model
std::string Model::tabularHeader(){
    std::string returnString = "Iteration\tPosterior\tLikelihood\tTree Prior";
    for(int i = 0; i < stateSpace; i++){
        returnString += "\tPi[" + std::to_string(i) + "]";
    }

    for(int i = 0; i < rateMatrix->transNameOrder().size(); i++){
        returnString += "\t" + rateMatrix->transNameOrder()[i];
    }

    return returnString + "\n";
}

// This method when called returns a string that prints out the values of each parameter. Parameters
// include each branch's gamma parameters, the stationary distribution, likelihood and etc.
std::string Model::tabularOut(int i){
    std::string returnString = std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" +
                               std::to_string(currentLikelihood) + "\t" + std::to_string(tree->lnPrior());
    std::vector<double> stationary = rateMatrix->getStationary();
    for(double i : stationary){
        returnString += "\t" + std::to_string(i);
    }

    // print nucleotide to nucleotide rates
    std::vector<double> rates = rateMatrix->getRate();
    for(int i = 0; i < rates.size(); i++){
        returnString += "\t" + std::to_string(rates[i]);
    }

    return returnString + "\n";
}

// call to get the header for tree output
std::string Model::treeHeader(){
    return "Iteration\tPosterior\tTree\n";
}

// call to get the tree output string for the current iteration for .tree
std::string Model::treeOut(int i){
    return std::to_string(i) + "\t" + std::to_string(lnPrior() + currentLikelihood) + "\t" + tree->writeNewick() + "\n";
}

