#include "TreeParameter.hpp"
#include "core/RandomVariable.hpp"
#include "TreeObject.hpp"
#include "core/Probability.hpp"
#include "core/Alignment.hpp"
#include "Node.hpp"
#include <cmath>
#include <test.h>

TreeParameter::TreeParameter(Alignment* aln) : currentPrior(0.0), oldPrior(0.0), 
                                                         moveChoice(-1), rateCount(0), rateAcceptCount(0), rateDelta(0.5), shapeCount(0),
                                                         shapeAcceptCount(0), shapeDelta(0.5), shapePriorRate(1.0), ratePriorRate(1.0), shapePriorShape(1.0), ratePriorShape(1.0){

    // first make a new tree based on aligment
    trees[0] = new TreeObject(aln);
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // iterate through each node and assign a branch gamma distribution (shape, rate), excluding the root
    std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
    for(Node* n : nodes) {
        if(n != trees[0]->getRoot()) {

            // shape and rate both have gamma priors
            double shape = Probability::Gamma::rv(&rng, shapePriorShape, shapePriorRate); 
            double rate = Probability::Gamma::rv(&rng, ratePriorShape, ratePriorRate); 
            trees[0]->setGammaDist(n, shape, rate);
        }
    }

    // copy the tree to the tree vector in case we reject
    trees[1] = new TreeObject(*trees[0]);

    // get the prior for all the branch gamma distributions
    double lnShape = 0.0;
    double lnRate = 0.0;
    std::vector<std::vector<double>> allGammaParams = trees[0]->getGammas();
    for(std::vector<double> gammaParam: allGammaParams){
        lnShape += Probability::Gamma::lnPdf(shapePriorShape, shapePriorRate, gammaParam[0]);
        lnRate += Probability::Gamma::lnPdf(shapePriorShape, ratePriorRate, gammaParam[1]);
    }
    currentPrior = lnShape + lnRate;
    dirty();
}

// delete the store tree objects
TreeParameter::~TreeParameter(){
    delete trees[0];
    delete trees[1]; 
}

// This method is in the event the proposed tree gets accepted, copy the accepted tree (tree[0]) to trees[1]
// for storage. Set currentPrior to oldPrior to store it. Also record what type of treeMove was made for the burn-in.
// Topology is untunable as its random.
void TreeParameter::accept(){
    *trees[1] = *trees[0]; // copy accepted tree to trees[1] for storage
    oldPrior = currentPrior;

    // update which acceptance occured
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
    *trees[0] = *trees[1]; // get the previous accepted tree back
    currentPrior = oldPrior;
    moveChoice = -1;
}

// This update randomly selects a node and then randomly selects either the shape or rate parameter of its
// ancestor branch to update.
double TreeParameter::updateTreeGamma(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

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
        double scale = std::exp(shapeDelta * (rng.uniformRv() - 0.5)); // propose new shape by scaling
        currentPrior -= Probability::Gamma::lnPdf(shapePriorShape, shapePriorRate, gammaParams[0]); // before rescaling, remove contribution to the prior of the current shape
        currentPrior += Probability::Gamma::lnPdf(shapePriorShape, shapePriorRate, gammaParams[0] * scale); // after rescaling, add contribution to the prior of the new shape               
        tree->setGammaDist(randNode, gammaParams[0] * scale, gammaParams[1]);
        hastings = std::log(scale);
    } else{
        rateCount += 1;
        moveChoice = 0;
        std::vector<double> gammaParams = tree->getGammaParams(randNode);
        double scale = std::exp(rateDelta * (rng.uniformRv() - 0.5));
        currentPrior -= Probability::Gamma::lnPdf(shapePriorShape, shapePriorRate, gammaParams[1]); // before rescaling, remove contribution to the prior of the current shape
        currentPrior += Probability::Gamma::lnPdf(shapePriorShape, shapePriorRate, gammaParams[1] * scale); // after rescaling, add contribution to the prior of the new shape               
        tree->setGammaDist(randNode, gammaParams[0], gammaParams[1] * scale);
        hastings = std::log(scale);
    }

    randNode->setNeedsTPUpdate(true); // branch gamma change requires a TP update
    this->dirty(); // mark tree as requiring update
    return hastings;
}

// This tree update changes the topology of the tree via NNI (nearest neighbor interchange)
// NNI implementation:
// pick a node p with the node a as its ancestor. Get p's "cousin" from node a and then swap the subtree with one of p's children
// Does not change anything related to branch lengths, only how the nodes are arranged
double TreeParameter::updateTreeMove() {
    RandomVariable& rng = RandomVariable::randomVariableInstance();
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
    while(internalNodeAncestor == nullptr || 
          internalNode == root || 
          internalNode->getIsTip()
    );

    // now we need to pick which of p's subtrees to swap
    std::set<Node*> iNNeighborSet = internalNode->getNeighbors();
    iNNeighborSet.erase(internalNodeAncestor);
    Node* swap1 = Node::chooseNodeFromSet(iNNeighborSet);

    // the other subtree to swap with is always p's cousin
    std::set<Node*> iNANeighborSet = internalNodeAncestor->getNeighbors();
    iNANeighborSet.erase(internalNode);
    iNANeighborSet.erase(internalNodeAncestor->getAncestor()); // this might not be safe
    Node* swap2 = Node::chooseNodeFromSet(iNANeighborSet); 

    // make swap1 a child of internalNodeAncestor by adding the relationship and removing its relationship with internalNode
    swap1->setAncestor(internalNodeAncestor);
    internalNodeAncestor->addNeighbor(swap1);
    swap1->addNeighbor(internalNodeAncestor);
    swap1->removeNeighbor(internalNode);
    internalNode->removeNeighbor(swap1);

    // make swap2 be a child of p
    swap2->setAncestor(internalNode);
    swap2->addNeighbor(internalNode);
    internalNode->addNeighbor(swap2);
    internalNodeAncestor->removeNeighbor(swap2);
    swap2->removeNeighbor(internalNodeAncestor);

    tree->initPostOrder(); // tree structure has changed, so reinitialize post-order node vector
    this->dirty(); // mark tree as requiring update for likelihood calculation
    double hastings = 0.0; // balanced move, so hastings == 0
    return hastings;
}

// During the burn-in phase of MCMC we seek to tune how the rate and shape parameter proposals are made for all
// nodes. The classic heurestic of aiming for an acceptance rate of ~0.33 is used. The topology itself can not be
// tuned, just the branch length param proposals (rate and shape).
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
    
    // Update shapeDelta the same way rateDelta is updated
    double shapeRate = (double)shapeAcceptCount/(double)shapeCount;
    if ( shapeRate > 0.33 ) {
        shapeDelta *= (1.0 + ((shapeRate-0.33)/0.67));
    }
    else {
        shapeDelta /= (2.0 - shapeRate/0.33);
    }

    // reset counts for the next round of tuning
    rateAcceptCount = 0;
    rateCount = 0;
    shapeAcceptCount = 0;
    shapeCount = 0;

    std::cout << "rateRate: " << rateRate << " shapeRate: " << shapeRate << "\n";
}
