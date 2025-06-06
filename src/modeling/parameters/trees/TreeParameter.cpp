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
                                                         treeCount(0), treeAcceptCount(0), treeAlpha(10000) {
    fixedTree = newick != ""; // fixedTree is true if newick is not empty and false otherwise   
    if(!fixedTree) // if newick is an empty string
        trees[0] = new TreeObject(aln);
    else // if newick is not empty this
        trees[0] = new TreeObject(newick, aln->getTaxaNames());

    // if not a fixed tree, meaning the topology
    if(!fixedTree){
        #ifdef TEST
        RandomVariable& rng = RandomVariable::randomVariableInstance(100);
        #endif
        #ifndef TEST
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        #endif
        std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
        // given that its postorder traversal, root_node will always be last
        for(Node* n : nodes) {
            // if the current node is the root
            if(n != trees[0]->getRoot()) {
                // a draw from the exponential distribution where it is a random rate and lambda is the time between occurences
                // idk what that really means I read that online and I am unsure why a draw from an exponential is the branch length
                // it may be that it is setting a flat prior for the branch length
                trees[0]->setBranchLength(n, Probability::Exponential::rv(&rng, lambda));
            }
        }
    }

    trees[1] = new TreeObject(*trees[0]);

    std::vector<double> values = trees[0]->getBranchLengths();
    double totalLength = 0.0;
    for(double val : values){
        totalLength += val;
    }
    currentPrior = Probability::Gamma::lnPdf(values.size(), lambda, totalLength);
    oldPrior = currentPrior;

    dirty();
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


double TreeParameter::update() {
    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(12);
    #endif
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    #endif
    double randomMove = rng.uniformRv();

    double hastings = 0.0;
    
    #ifndef TEST
    if(randomMove < 0.75){
    #endif
     #ifdef TEST
    if(randomMove > 0){
    #endif
        
        // Change topology and branch lengths because it is not a fixedTree
        if(!fixedTree){

            // NNI implementation
            // picks a branch containing subtrees s1, s2, s3 and s4 in the configuration ((s1, s2), s3, s4)
            // and randomly transforms the branch into either ((s1, s3), s2, s4) or ((s1, s4), s2, s3). Swapping
            // out an internal subtree with a subtree that diverged earlier
            if(1){

                #ifdef TEST
                std::cout << "moving NNI\n";
                #endif
                moveChoice = 2; //2 represents NNI cause I said so
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
                while(internalNodeAncestor == nullptr || internalNodeAncestor == root || internalNode->getIsTip());

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
                    if(internalNodeAncestorAncestor == root && tempnode != internalNodeAncestor){
                        internalNodeAncestorSibling = tempnode;
                        break;
                    }
                    else if (tempnode != internalNodeAncestorAncestor && tempnode != internalNodeAncestor ){
                        internalNodeAncestorSibling = tempnode;
                        break;
                    }
                }
                Node *s4 = internalNodeAncestorSibling;

                // now that we have our 4 subtrees, we need to do a coin flip to decide which rearrangement
                // to use. Fundamentally the swaps are the same but its just with different subtrees
                int coinFlip = (int)(rng.uniformRv() * 1);
                Node *swap1 = s2;
                Node *swap2;
                
                // swap1 = s2 and swap2 = s3
                if(coinFlip){
                    swap2 = s3;
                } 

                // swap1 = s2 and swap2 = s4
                else{
                    swap2 = s4;
                }
                // we are essentially swapping swap1 and swap2, store swap1's neighbor and ancestor
                // before setting swap1's to swap 2, then swap 2 to swap 1
                Node* swap1Ancestor = swap1->getAncestor();
                std::set<Node*> swap1NeighborSet = swap1->getNeighbors();
                std::set<Node*> swap2NeighborSet = swap2->getNeighbors();
                std::set<Node*> storeset;
                for(Node* n : swap1NeighborSet){
                    storeset.insert(n);
                }
                swap1->removeAllNeighbors();

                // change swap1's neighbors and ancestors to swap2
                swap1->setAncestor(swap2->getAncestor());
                for(Node* n : swap2NeighborSet){
                    swap1->addNeighbor(n);
                }

                // now "move" swap2 into swap1's location
                swap2->removeAllNeighbors();
                swap2->setAncestor(swap1Ancestor);
                for(Node* n: storeset){
                    swap2->addNeighbor(n);
                }

                // now go through swap1's ancestors and neighbors and remove any mention of swap2
                // and replace it with itself
                for(Node *n: swap2NeighborSet){
                    n->removeNeighbor(swap2);
                    n->addNeighbor(swap1);
                    Node *ancestor = n->getAncestor();
                    if(ancestor == swap2){
                        n->setAncestor(swap1);
                    }
                }

                // now do the same for swap2
                for(Node *n: swap1NeighborSet){
                    n->removeNeighbor(swap1);
                    n->addNeighbor(swap2);
                    Node *ancestor = n->getAncestor();
                    if(ancestor == swap1){
                        n->setAncestor(swap2);
                    }
                }
                                
                    
                
            } 
            else{
                moveChoice = 0;
                branchCount += 1;
                TreeObject* tree = trees[0];
                std::vector<Node*> nodes = tree->getPostOrderSeq();
                Node* root = tree->getRoot();

                Node* u = nullptr;
                do{
                    u = nodes[(int)(rng.uniformRv() * nodes.size())];
                }
                while(u == root || u->getIsTip() == true);
                Node* v = u->getAncestor();

                std::set<Node*> neighbors1 = u->getNeighbors();
                neighbors1.erase(v);//Exclude v
                Node* a = Node::chooseNodeFromSet(neighbors1);

                std::set<Node*> neighbors2 = v->getNeighbors();
                neighbors2.erase(u);//Don't select u
                Node* c = Node::chooseNodeFromSet(neighbors2);

                double scale = std::exp(branchDelta * (rng.uniformRv() - 0.5));

                double paths[3];
                paths[0] = tree->getBranchLength(u) * scale;
                paths[1] = tree->getBranchLength(a) * scale;
                Node* b3 = nullptr;
                if(c != v->getAncestor())
                    b3 = c;
                else
                    b3 = v;
                paths[2] = tree->getBranchLength(b3) * scale;

                double totalPath = paths[0] + paths[1] + paths[2];

                std::vector<Node*> nodeSet = {a, b3};
                int pick = (int)(rng.uniformRv() * 2);
                double randomLoc = rng.uniformRv() * totalPath;

                //The pick decides the oritentation of the path
                if(randomLoc <= totalPath - paths[pick]){
                    tree->setBranchLength(nodeSet[pick], randomLoc);
                    tree->setBranchLength(u, totalPath - paths[pick] - randomLoc);
                    tree->setBranchLength(nodeSet[pick ^ 1], paths[pick]);
                }
                else{
                    u->removeNeighbor(a);
                    a->removeNeighbor(u);
                    v->removeNeighbor(c);
                    c->removeNeighbor(v);

                    v->addNeighbor(a);
                    a->addNeighbor(v);
                    u->addNeighbor(c);
                    c->addNeighbor(u);
                    tree->setBranchLength(nodeSet[pick ^ 1], totalPath - randomLoc);
                    tree->setBranchLength(u, randomLoc - (totalPath - paths[pick]));
                    tree->setBranchLength(nodeSet[pick], totalPath - paths[pick]);

                    //Rooting logic
                    if(v->getAncestor() == c){
                        u->setAncestor(c);
                        v->setAncestor(u);
                        a->setAncestor(v);
                    }
                    else{
                        c->setAncestor(u);
                        u->setAncestor(v);
                        a->setAncestor(v);
                    }
                }


                u->setNeedsTPUpdate(true);
                v->setNeedsTPUpdate(true);
                a->setNeedsTPUpdate(true);
                c->setNeedsTPUpdate(true);

                //The updating gets a little awkward because we don't really know the branching here.
                Node* q = v;
                if(v->getAncestor() != u)
                    q = u;

                do{
                    if(q->getIsTip() == false)
                        q->setNeedsCLUpdate(true);
                    q = q->getAncestor();
                }
                while(q != root);
                root->setNeedsCLUpdate(true);

                tree->initPostOrder();
                this->dirty();

                hastings = 3 * std::log(scale);
                }

            }

            else {
                moveChoice = 0;
                branchCount += 1;
                std::vector<Node*> nodes = trees[0]->getPostOrderSeq();
                Node* root = trees[0]->getRoot();

                Node* p = nullptr;
                do{
                    p = nodes[(int)(rng.uniformRv() * nodes.size())];
                }
                while(p == root);

                double currentV = trees[0]->getBranchLength(p);
                double scale = std::exp(branchDelta * (rng.uniformRv() - 0.5));
                double newV = currentV * scale;
                trees[0]->setBranchLength(p, newV);
                p->setNeedsTPUpdate(true);

                Node* q = p;
                do{
                    if(q->getIsTip() == false)
                        q->setNeedsCLUpdate(true);
                    
                    q = q->getAncestor();
                } 
                while(q != root);
                root->setNeedsCLUpdate(true);

                this->dirty();

                hastings = std::log(scale);
        }
    }
    else {
        moveChoice = 1;
        treeCount += 1;
        std::map<Node*, double> branchMapping = trees[0]->getBranchLengthMapping();
        trees[0]->updateAll();
        this->dirty();

        std::vector<double> values;
        std::vector<Node*> nodeIndices;
        double totalLength = 0.0;
        for(auto mapping : branchMapping){
            double l = mapping.second;
            nodeIndices.push_back(mapping.first);
            values.push_back(l);
            totalLength += l;
        }

        std::vector<double> alphaForward(values.size(), 0.0);
        std::vector<double> alphaReverse(values.size(), 0.0);
        std::vector<double> z(values.size(), 0.0);

        for(int i = 0; i < values.size(); i++) {
            values[i] /= totalLength;
            alphaForward[i] = values[i] * treeAlpha;
        }
        
        Probability::Dirichlet::rv(&rng, alphaForward, z);

        for(int i = 0; i < z.size(); i++) {
            alphaReverse[i] = z[i] * treeAlpha;
        }
        
        hastings  = Probability::Dirichlet::lnPdf(alphaReverse, values) - Probability::Dirichlet::lnPdf(alphaForward, z);

        for(int i = 0; i < values.size(); i++){
            trees[0]->setBranchLength(nodeIndices[i], z[i] * totalLength);
        } 
    }

    std::vector<double> values = trees[0]->getBranchLengths();
    double totalLength = 0.0;
    for(double val : values){
        totalLength += val;
    }
    currentPrior = Probability::Gamma::lnPdf(values.size(), lambda, totalLength);

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