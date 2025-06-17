#include "TreeParameter.hpp"
#include "core/RandomVariable.hpp"
#include "TreeObject.hpp"
#include "core/Probability.hpp"
#include "core/Alignment.hpp"
#include "Node.hpp"
#include <cmath>
#include <test.h>

TreeParameter::TreeParameter(Alignment* aln, std::string newick, double l) : lambda(l), currentPrior(0.0), oldPrior(0.0), 
                                                         branchDelta(1), moveChoice(-1), branchCount(0), branchAcceptCount(0), 
                                                         treeCount(0), treeAcceptCount(0), treeAlpha(10000), gammaDelta(1.1) {
    fixedTree = newick != ""; // fixedTree is true if newick is not empty and false otherwise   
    if(!fixedTree) // if newick is an empty string
        trees[0] = new TreeObject(aln);
    // else // if newick is not empty this
    //     trees[0] = new TreeObject(newick, aln->getTaxaNames());

    // if not a fixed tree, meaning the topology
    if(!fixedTree){
        #ifdef TEST
        RandomVariable& rng = RandomVariable::randomVariableInstance(100);
        #endif
        #ifndef TEST
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        #endif
        std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
        for(Node* n : nodes) {
            if(n != trees[0]->getRoot()) {
                double alpha = Probability::Gamma::rv(&rng, 1, 1);
                double beta = Probability::Gamma::rv(&rng, 1, 1);
                trees[0]->setGammaDist(n, 1, 1);
            }
        }
    }

    trees[1] = new TreeObject(*trees[0]);

    // FIX THIS BY IMPLEMENTING THE NEW GAMMA LIKELIHOOD CALCULATION
    // std::vector<double> values = trees[0]->getBranchLengths();
    // double totalLength = 0.0;
    // for(double val : values){
    //     totalLength += val;
    // }
    // currentPrior = Probability::Gamma::lnPdf(values.size(), lambda, totalLength);
    // oldPrior = currentPrior;

    dirty();

    #ifdef TEST
    trees[0]->setNodeNameIndex();
    #endif
}

TreeParameter::~TreeParameter(){
    delete trees[0];
    delete trees[1]; 
}

void TreeParameter::accept(){
    *trees[1] = *trees[0];
    oldPrior = currentPrior;

    if(moveChoice == 0){
        branchAcceptCount += 1;
    }
    else if(moveChoice == 1){
        treeAcceptCount += 1;
    }

    moveChoice = -1;
}

void TreeParameter::reject(){
    *trees[0] = *trees[1];
    currentPrior = oldPrior;
    moveChoice = -1;
}

// IMPLEMENT THE ABILITY TO STEP AND A PROPER RESCALING 
double TreeParameter::updateTreeGamma(){
    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(12);
    #endif
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    #endif
    double hastings = 0.0;

    moveChoice = 1;
    treeCount += 1;

    // pick a random node that is not the root to update one of their ancestor branch's gamma parameters
    TreeObject* tree = trees[0];
    std::vector<Node*> nodes = tree->getPostOrderSeq();
    Node* randNode = nullptr;
    do{
        randNode = nodes[(int)(rng.uniformRv() * nodes.size())];
    } while(randNode != tree->getRoot());
    
    // flip a coin to determine which parameter to update (alpha if 1 or beta if 0)
    int coinFlip = (int)(rng.uniformRv() * 1);
    std::vector<double> gammaParams = tree->getGammaParams(randNode);
    double changeParam = gammaParams[coinFlip];
    
    // now do a simple rescaling proposal by drawing a multiplicative factor from a normal distribution
    // CHECK TO SEE WHAT THE MU AND SIGMA OF NORMAL RV is
    // double scale = std::exp(branchDelta * (rng.uniformRv() - 0.5));
    double scale = Probability::Normal::rv(&rng);
    double proposeParam = changeParam * gammaDelta;

    // the hastings is assymetrical as it is multiplicative, so lnpdf (reverse move) - lnpdf (forward move)
    hastings = Probability::Normal::lnPdf(0, 1 , 1/gammaDelta) - Probability::Normal::lnPdf(0, 1, gammaDelta);


    // std::map<Node*, std::vector<double>> gammaMapping = trees[0]->getGammaMap();
    // trees[0]->updateAll();
    // this->dirty();

    // std::vector<double> changeParam;
    // std::vector<double> otherParam;
    // std::vector<Node*> nodeIndices;
    // double totalSum = 0.0;
    
    // // there are two gammaParams to do an update for, pick one randomly
    // // we will propose an update to alpha (shape) if coinFlip is true or to beta (shape) is false
    // int coinFlip = (int)(rng.uniformRv() * 1);
    // int changeParamIndex = coinFlip ? 0 : 1;
    // int otherParamIndex = coinFlip ? 1 : 0;
    // for(auto mapping : gammaMapping){
    //     double change = mapping.second[changeParamIndex];
    //     double other = mapping.second[otherParamIndex];
    //     nodeIndices.push_back(mapping.first);
    //     changeParam.push_back(change);
    //     otherParam.push_back(other);
    //     totalSum += change;
    // }

    // // create some empty vectors filled with 0
    // std::vector<double> alphaForward(changeParam.size(), 0.0);
    // std::vector<double> alphaReverse(changeParam.size(), 0.0);
    // std::vector<double> z(changeParam.size(), 0.0);

    // // sum normalize the dataset
    // for(int i = 0; i < changeParam.size(); i++) {
    //     changeParam[i] /= totalSum;
    //     alphaForward[i] = changeParam[i] * treeAlpha; //NOTE: may need to change treeAlpha later so it mixes better
    // }
    
    // Probability::Dirichlet::rv(&rng, alphaForward, z);

    // for(int i = 0; i < z.size(); i++) {
    //     alphaReverse[i] = z[i] * treeAlpha; // how come we don't divide by treeAlpha here
    // }
    
    // hastings  = Probability::Dirichlet::lnPdf(alphaReverse, changeParam) - Probability::Dirichlet::lnPdf(alphaForward, z);

    // for(int i = 0; i < changeParam.size(); i++){
    //     if(coinFlip){
    //         trees[0]->setGammaDist(nodeIndices[i], z[i]*totalSum, otherParam[i]);
    //     } else{
    //         trees[0]->setGammaDist(nodeIndices[i], otherParam[i], z[i]*totalSum);
    //     }
    // } 

}

double TreeParameter::updateTreeMove() {
    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(12);
    #endif
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    #endif
    double hastings = 0.0;

    // NNI implementation
    // picks a branch containing subtrees s1, s2, s3 and s4 in the configuration ((s1, s2), s3, s4)
    // and randomly transforms the branch into either ((s1, s3), s2, s4) or ((s1, s4), s2, s3). Swapping
    // out an internal subtree with a subtree that diverged earlier
    moveChoice = 0; 
    branchCount += 0;
    TreeObject* tree = trees[0];
    std::vector<Node*> nodes = tree->getPostOrderSeq();
    Node* root = tree->getRoot();

    // pick a random internal branch by picking a node's whose ancestor branch
    // is an internal branch. The additional stipulation is that the ancestor of the internalNode
    // (internalNodeAncestor) cannot be the root
    Node* internalNode = nullptr;
    Node* internalNodeAncestor = nullptr;
    do{
        internalNode = nodes[(int)(rng.uniformRv() * nodes.size())];

        // if a node isn't the root, it will always have an ancestor
        if(internalNode != root){
            internalNodeAncestor = internalNode->getAncestor();
        }
    }
    while(internalNodeAncestor == nullptr || internalNodeAncestor == root || internalNode == root || internalNode->getIsTip());

    // we know that internalNodeAncestor is NOT the root, so we can safely get its ancestor
    Node* internalNodeAncestorAncestor = internalNodeAncestor->getAncestor();

    // the descendants of the internalNode are both subtrees s1 and s2
    std::set<Node*> internalNodeSet = internalNode->getNeighbors();
    Node* s1 = nullptr;
    Node* s2 = nullptr;
    Node* tempnode = nullptr;

    // keep picking nodes from the internalNode's neighbors (internalNodeSet)
    // until we populate s1 and s2 with the internalNode's descendants
    while(1){
        tempnode = tempnode->chooseNodeFromSet(internalNodeSet);
        if(tempnode != internalNodeAncestor){
            if(s1 == nullptr){
                s1 = tempnode;
            } 
            else if(s1 != tempnode && s2 == nullptr){
                s2 = tempnode;
                break;
            }
        }
    }

    // internalNode has a sibling (the other descendant of internalNodeAncestor)
    // get internalNode's sibling as it is s3
    Node* internalSiblingNode = nullptr;
    std::set<Node*> internalNodeAncestorNeighborSet = internalNodeAncestor->getNeighbors();
    while(1){
        tempnode = tempnode->chooseNodeFromSet(internalNodeAncestorNeighborSet);
        if(tempnode != internalNodeAncestorAncestor && tempnode != internalNode){
            internalSiblingNode = tempnode;
            break;
        }
    }
    Node *s3 = internalSiblingNode;

    // s4 is the sibling of internalNodeAncestor, so repeat same process again but take
    // into consideration that internalNodeAncestorAncestor may be the root
    Node * internalNodeAncestorSibling = nullptr;
    std::set<Node*> internalNodeAncestorAncestorNeighborSet = internalNodeAncestorAncestor->getNeighbors();
    while(1){
        tempnode = tempnode->chooseNodeFromSet(internalNodeAncestorAncestorNeighborSet);
        if (tempnode != internalNodeAncestorAncestor->getAncestor() && tempnode != internalNodeAncestor ){
            internalNodeAncestorSibling = tempnode;
            break;
        }
    }
    Node *s4 = internalNodeAncestorSibling;

    // now that we have our 4 subtrees, we need to do a coin flip to decide which rearrangement
    // to use. Fundamentally the swaps are the same but its just with different subtrees
    int coinFlip = (int)(rng.uniformRv() * 1);
    Node *swap1 = s2;
    Node *swap2 = coinFlip ? s3 : s4;

    #ifdef TEST
        std::cout << "\nswapping " << swap1->getIndex() << " and " << swap2->getIndex() << "\n" << std::flush;
    #endif

    Node* swap1Ancestor = swap1->getAncestor();
    Node* swap2Ancestor = swap2->getAncestor();
    swap1->removeNeighbor(swap1Ancestor);
    swap1->addNeighbor(swap2Ancestor);
    swap1->setAncestor(swap2Ancestor);
    swap2->removeNeighbor(swap2Ancestor);
    swap2->addNeighbor(swap1Ancestor);
    swap2->setAncestor(swap1Ancestor);

    // now go to the ancestor's of swap1 and swap2 and make them point to their new children
    swap1Ancestor->removeNeighbor(swap1);
    swap1Ancestor->addNeighbor(swap2);
    swap2Ancestor->removeNeighbor(swap2);
    swap2Ancestor->addNeighbor(swap1);

    // set some flags for the nodes affected by the changes, basically all nodes of the subtrees that 
    // got swapped and the flow via ancestors back to the root need to have CL update
    Node *needsCLupdate = swap1;
    while(needsCLupdate != root){
        if(!needsCLupdate->getIsTip()){
            needsCLupdate->setNeedsCLUpdate(true);
        }
        needsCLupdate = needsCLupdate->getAncestor();
    }

    needsCLupdate = swap2;
    while(needsCLupdate != root){
        if(!needsCLupdate->getIsTip()){
            needsCLupdate->setNeedsCLUpdate(true);
        }
        needsCLupdate = needsCLupdate->getAncestor();
    }
    root->setNeedsCLUpdate(true);

    // tree flags and hastings time
    tree->initPostOrder();
    this->dirty();
    hastings = 0; // we are equally likely to go back to where we started intuitively

    return hastings;
}

void TreeParameter::tune() {
    double rate1 = (double)branchAcceptCount/(double)branchCount;

    if ( rate1 > 0.33 ) {
        branchDelta *= (1.0 + ((rate1-0.33)/0.67));
    }
    else {
        branchDelta /= (2.0 - rate1/0.33);
    }
    branchAcceptCount = 0;
    branchCount = 0;

    double rate2 = (double)treeAcceptCount/(double)treeCount;

    if ( rate2 > 0.33 ) {
        treeAlpha /= (1.0 + ((rate2-0.33)/0.67));
    }
    else {
        treeAlpha *= (2.0 - rate2/0.33);
    }
    treeAcceptCount = 0;
    treeCount = 0;
}

double TreeParameter::lnPrior() {
    return currentPrior;
}