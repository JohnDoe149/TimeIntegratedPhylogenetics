#include "TreeParameter.hpp"
#include "core/RandomVariable.hpp"
#include "TreeObject.hpp"
#include "core/Probability.hpp"
#include "core/Alignment.hpp"
#include "Node.hpp"
#include <cmath>
#include <test.h>

TreeParameter::TreeParameter(Alignment* aln, std::string newick, double l) : lambda(l), currentPrior(0.0), oldPrior(0.0), 
                                                         moveChoice(-1), rateCount(0), rateAcceptCount(0), rateDelta(0.5), shapeCount(0),
                                                         shapeAcceptCount(0), shapeDelta(0.5){

    // first make a new tree based on aln
    trees[0] = new TreeObject(aln);
    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(100);
    #endif
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    #endif

    // Now go through each node and assign their ancestor branch random gamma parameters shape and rate.
    // except the root as it has no ancestor branch
    std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
    for(Node* n : nodes) {
        if(n != trees[0]->getRoot()) {

            // alpha and beta both have gamma priors
            double shape = Probability::Gamma::rv(&rng, 4, 10);
            double rate = Probability::Gamma::rv(&rng, 4, 10);
            trees[0]->setGammaDist(n, shape, rate);
        }
    }

    // copy the tree to the tree vector in case we reject
    trees[1] = new TreeObject(*trees[0]);

    // go through each node and calculate the probability of the gamma parameters and add them up 
    // to use in the posterior calculation
    double lnShape = 0.0;
    double lnRate = 0.0;
    std::vector<std::vector<double>> allGammaParams = trees[0]->getGammas();
    for(std::vector<double> gammaParam: allGammaParams){
        lnShape += Probability::Gamma::lnPdf(4, 10, gammaParam[0]);
        lnRate += Probability::Gamma::lnPdf(4, 10, gammaParam[1]);
    }
    currentPrior = lnShape + lnRate;
    dirty();

    #ifdef TEST
    trees[0]->setNodeNameIndex();
    #endif
}

TreeParameter::~TreeParameter(){
    delete trees[0];
    delete trees[1]; 
}

// This method is in the event the proposed tree gets accepted, copy the accepted tree (tree[0]) to trees[1]
// for storage. Set currentPrior to oldPrior to store it. Also record what type of treeMove was made for the burn-in
void TreeParameter::accept(){
    *trees[1] = *trees[0];
    oldPrior = currentPrior;

    if(moveChoice == 0){
        rateAcceptCount += 1;
    } 
    else if(moveChoice == 1){
        shapeAcceptCount += 1;
    }

    // reset moveChoice
    moveChoice = -1;
}

// In the event the proposed tree is rejected, copy the tree we stored in trees[1] to use for the next proposal
void TreeParameter::reject(){
    *trees[0] = *trees[1];
    currentPrior = oldPrior;
    moveChoice = -1;
}

// This update randomly selects a node and then randomly selects either the shape or rate parameter of its
// ancestor branch to update.
double TreeParameter::updateTreeGamma(){
    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(12);
    #endif
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    #endif

    // pick a random node that is not the root to update one of their ancestor branch's gamma parameters
    TreeObject* tree = trees[0];
    std::vector<Node*> nodes = tree->getPostOrderSeq();
    Node* randNode = nullptr;
    Node* treeRoot = tree->getRoot();
    do{
        randNode = nodes[(int)(rng.uniformRv() * nodes.size())];
    } while(randNode == treeRoot);
    
    // flip a coin to determine which parameter to update (shape if 1 or rate if 0)
    int coinFlip = rng.uniformRv() < 0.5 ? 0 : 1;
    double hastings;

    // record whether or not a shape or rate update is going to be made
    if(coinFlip){
        shapeCount += 1;
        moveChoice = 1;
        std::vector<double> gammaParams = tree->getGammaParams(randNode);
        double scale = std::exp(shapeDelta * (rng.uniformRv() - 0.5));
        tree->setGammaDist(randNode, gammaParams[0] * scale, gammaParams[1]);
        hastings = std::log(scale);
    } else{
        rateCount += 1;
        moveChoice = 0;
        std::vector<double> gammaParams = tree->getGammaParams(randNode);
        double scale = std::exp(rateDelta * (rng.uniformRv() - 0.5));
        tree->setGammaDist(randNode, gammaParams[0], gammaParams[1] * scale);
        hastings = std::log(scale);
    }

    // set flags for the changed node, a changed shape or rate means the transition probabilty changes and all nodes
    // it "descended" from back to the root need to have their conditional likelihood recalculated
    // due to how Felsenstein's algorithm is calculated (from descendant to ancestor)
    randNode->setNeedsTPUpdate(true);
    if(!randNode->getIsTip()){
        randNode->setNeedsCLUpdate(true);
    }
        Node* randNodeAnc = randNode->getAncestor();
    while(randNodeAnc != treeRoot){
        randNodeAnc->setNeedsCLUpdate(true);
        randNodeAnc = randNodeAnc->getAncestor();
    }
    treeRoot->setNeedsCLUpdate(true);
    tree->initPostOrder();
    this->dirty();

    // go through each node and calculate the probability of the gamma parameters and add them up
    // to get our prior probability for the posterior calculation
    double lnShape = 0.0;
    double lnRate = 0.0;
    std::vector<std::vector<double>> allGammaParams = trees[0]->getGammas();
    for(std::vector<double> gammaParam: allGammaParams){
        lnShape += Probability::Gamma::lnPdf(4, 10, gammaParam[0]);
        lnRate += Probability::Gamma::lnPdf(4, 10, gammaParam[1]);
    }
    currentPrior = lnShape + lnRate;

    return hastings;
}

// This tree update changes the topology of the tree via NNI (nearest neighbor interchange)
// NNI implementation:
// picks a branch containing subtrees s1, s2, s3 and s4 in the configuration ((s1, s2), s3, s4)
// and randomly transforms the branch into either ((s1, s3), s2, s4) or ((s1, s4), s2, s3). Swapping
// out an internal subtree with a subtree that diverged earlier
// Does not change anything related to branch lengths, only how the nodes are arranged
double TreeParameter::updateTreeMove() {
    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(12);
    #endif
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    #endif

    double hastings;
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
    int coinFlip = rng.uniformRv() < 0.5 ? 0 : 1;
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

// During the burn-in phase of MCMC we seek to tune how the rate and shape parameter proposals are made for all
// nodes. The classic heurestic of aiming for an acceptance rate of ~0.33 is used. The topology itself can not be
// tuned, just the branch lengths (rate and shape).
void TreeParameter::tune() {

    // Calculate the rate proposals acceptance rate, if its too high, make rateDelta larger to decrease acceptance
    // otherwise make it smaller (more conservative proposals) to increase acceptance.
    double rateRate = (double)rateAcceptCount/(double)rateCount;
    if ( rateRate > 0.33 ) {
        rateDelta *= (1.0 + ((rateRate-0.33)/0.67));
    }
    else {
        rateDelta /= (2.0 - rateRate/0.33);
    }
    rateAcceptCount = 0;
    rateCount = 0;
    
    // Update shapeDelta the same way rateDelta is updated
    double shapeRate = (double)shapeAcceptCount/(double)shapeCount;
    if ( shapeRate > 0.33 ) {
        shapeDelta *= (1.0 + ((shapeRate-0.33)/0.67));
    }
    else {
        shapeDelta /= (2.0 - shapeRate/0.33);
    }
    shapeAcceptCount = 0;
    shapeCount = 0;

    std::cout << "rateRate: " << rateRate << " shapeRate: " << shapeRate << "\n";
}

// A simple "getter" method to get the prior for a tree proposal 
double TreeParameter::lnPrior() {
    return currentPrior;
}