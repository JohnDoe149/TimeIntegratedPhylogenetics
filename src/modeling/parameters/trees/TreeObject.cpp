#include "TreeObject.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "Node.hpp"
#include "core/Msg.hpp"
#include "core/Alignment.hpp"
#include <iostream>
#include <cmath>
#include "test.h"

/*
=======================================================================
                    TREE CONSTRUCTORS/DESTRUCTORS
=======================================================================
*/

TreeObject::TreeObject(int nt) : numTaxa(nt) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    root = addNode();
    root->setName("Root");
    for(int i = 0; i < 2; i++) {
        Node* p = addNode();
        p->setName("Taxon_" + std::to_string(i+1));
        p->setIsTip(true);
        p->setIndex(i);
        p->addNeighbor(root);
        root->addNeighbor(p);
        p->setAncestor(root);
    }

    // build up the full tree by randomly adding branches to the existing tree
    for(int i = 2; i < numTaxa; i++){
        Node* p = nullptr;
        do {
            p = nodes[(int)(rng.uniformRv() * nodes.size())];
        } while(p == root);

        Node* pAnc = p->getAncestor();
        Node* newTip = addNode();
        newTip->setIsTip(true);
        newTip->setIndex(i);
        newTip->setName("Taxon_" + std::to_string(i+1));
        Node* newAnc = addNode();

        p->removeNeighbor(pAnc);
        p->addNeighbor(newAnc);
        p->setAncestor(newAnc);

        newAnc->addNeighbor(p);
        newAnc->addNeighbor(newTip);
        newAnc->addNeighbor(pAnc);
        newAnc->setAncestor(pAnc);

        newTip->addNeighbor(newAnc);
        newTip->setAncestor(newAnc);

        pAnc->removeNeighbor(p);
        pAnc->addNeighbor(newAnc);
    }
        
    // initialize the down pass sequence
    initPostOrder();

    // index the interior nodes (the tip nodes are indexed, above)
    int intIdx = numTaxa;
    for (int i=0, n=(int)postOrderSeq.size(); i<n; i++) {
        Node* p = postOrderSeq[i];
        if (p->getIsTip() == false)
            p->setIndex(intIdx++);
    }
    std::vector<double> zero_vec(numTaxa*2-2, 1.0);
    branchLengths = zero_vec;
}

//Generate a random tree and connect it to an alignment
TreeObject::TreeObject(Alignment* aln) : TreeObject(aln->getNumTaxa()) {

    std::vector<std::string> names = aln->getTaxaNames();
    for(Node* n : postOrderSeq){
        if(n->getIsTip()){
            n->setName(names[n->getIndex()]);
        }
    }
}

TreeObject::TreeObject(const TreeObject& t){
    clone(t);
}

TreeObject::~TreeObject(void) {
    deleteAllNodes();
}

//Deep Copy Operation
TreeObject& TreeObject::operator=(const TreeObject& rhs){
    if(this == &rhs)
        return *this;
    
    clone(rhs);
    return *this;
}

/*
=======================================================================
                        Other Functions
=======================================================================
*/

Node* TreeObject::addNode(void) {

    Node* newNode = new Node;
    newNode->setOffset((int)nodes.size());
    nodes.push_back(newNode);
    return newNode;
}

void TreeObject::clone(const TreeObject& t){
    if(nodes.size() != t.nodes.size()){
        deleteAllNodes();
        for(int i = 0; i < t.nodes.size(); i++)
            addNode();
    }
    this->branchGamma.clear();

    this->numTaxa = t.numTaxa;
    this->root = this->nodes[t.root->getOffset()];

    for(int i = 0; i < t.nodes.size(); i++){
        Node* p = this->nodes[i];
        Node* q = t.nodes[i];
        p->setIndex(q->getIndex());
        p->setIsTip(q->getIsTip());
        p->setName(q->getName());
        p->setNeedsCLUpdate(q->getNeedsCLUpdate());
        p->setNeedsTPUpdate(q->getNeedsTPUpdate());

        p->removeAllNeighbors();
        std::set<Node*>& qNeighbors = q->getNeighbors();
        for(Node* n : qNeighbors)
            p->addNeighbor(this->nodes[n->getOffset()]);

        if(q->getAncestor() != nullptr){
            Node* ancestor = this->nodes[q->getAncestor()->getOffset()];
            std::vector<double> gp = t.getGammaParams(q);
            p->setAncestor(ancestor);
            this->setGammaDist(p, gp[0], gp[1]);
        }
        else
            p->setAncestor(nullptr);
    }

    this->postOrderSeq.clear();
    for(int i = 0; i < t.postOrderSeq.size(); i++){
        Node* p = t.postOrderSeq[i];
        this->postOrderSeq.push_back(this->nodes[p->getOffset()]);
    }
}

void TreeObject::deleteAllNodes(){
    for (int i = 0; i < nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

void TreeObject::setGammaDist(Node* n, double shape, double rate){
    std::vector<double> gammaParams;
    gammaParams.push_back(shape);
    gammaParams.push_back(rate);
    auto it = branchGamma.find(n);

    if(it == branchGamma.end())
        branchGamma.insert(std::make_pair(n, gammaParams));
    else
        it->second = gammaParams;
}


std::vector<double> TreeObject::getGammaParams(Node* n) const{
    auto it = branchGamma.find(n);

    if(it == branchGamma.end())
        Msg::error("Couldn't find GammaParams for this node");
    return it->second;
}

std::vector<std::vector<double>> TreeObject::getGammas(){
    std::vector<std::vector<double>> returnVec;
    returnVec.reserve(branchGamma.size());

    for (auto v: branchGamma){

        // load all the gamma vectors into returnVec
        returnVec.push_back(v.second);
    }
    return returnVec;
}

std::map<Node*, std::vector<double>> TreeObject::getGammaMap(){
    return branchGamma;
}

std::string TreeObject::getNewick() const{
    std::stringstream strm;
    writeNode(root, strm);
    return strm.str();
}

std::vector<Node*> TreeObject::getTips() {
    std::vector<Node*> out;
    out.reserve(numTaxa);

    for(Node* n : nodes){
        if(n->getIsTip())
            out.push_back(n);
    }

    return out;
}

int TreeObject::getTaxonIndex(std::string token, std::vector<std::string> taxaNames){

    for(int i = 0, n = taxaNames.size(); i < n; i++){
        if(taxaNames[i] == token)
            return i;
    }

    return -1;
}

void TreeObject::initPostOrder(void) {
    postOrderSeq.clear();
    passDown(root, postOrderSeq);
}

std::vector<std::string> TreeObject::parseNewickString(std::string newick){
    std::vector<std::string> tokens;
    std::string str = "";
    for(int i = 0; i < newick.length(); i++){
        char c = newick[i];
        if(c == '(' || c == ')' || c == ',' || c == ':' || c == ';'){
            if(str != ""){
                tokens.push_back(str);
                str = "";
            }

            tokens.push_back(std::string(1, c));
        }
        else {
            str += std::string(1, c);
        }
    }

    return tokens;
}

void TreeObject::passDown(Node* p, std::vector<Node*>& vec) {

    if(p == nullptr)
        return;
    
    std::set<Node*>& pNeighbors = p->getNeighbors();

    for(Node* n : pNeighbors){
        if(n != p->getAncestor())
            passDown(n, vec);
    }

    vec.push_back(p);
}

void TreeObject::print(std::string header) const{
    std::cout << header << std::endl;
    print();
}

void TreeObject::print(void) const{

    showNode(root, 0);
}

//Output nodes of the tree in a whitespace-indented format
void TreeObject::showNode(Node* p, int indent) const{

    if(p == nullptr)
        return;

    for(int i = 0; i < indent; i++)
        std::cout << " ";

    std::cout << p->getIndex() << " ( ";
    std::set<Node*>& pNeighbors = p->getNeighbors();
    for (Node* d : pNeighbors)
        {
        if (d == p->getAncestor())
            std::cout << "a_";
        std::cout << d->getIndex() << " ";
        }
    std::cout << ") ";

    if(p->getAncestor() != nullptr)
        // FIX THIS LATER
        //<< this->getBranchLength(p) << " ";

    std::cout << p->getName();

    if (p == root)
        std::cout << " <- Root";
    std::cout << std::endl;

    for(Node* n : pNeighbors)
        {
        if(n != p->getAncestor())
            showNode(n, indent+3);
        }
}

void TreeObject::updateAll(){
    for(Node* n : nodes){
        if(n->getIsTip() == false)
            n->setNeedsCLUpdate(true);
        n->setNeedsTPUpdate(true);
    }
}

//special helper method for making all node's name their index for easy debugging
void TreeObject::setNodeNameIndex(){
    initPostOrder();

    // Initialize branch lengths
    for (int i=0, n=(int)postOrderSeq.size(); i<n; i++) {
        Node* p = postOrderSeq[i];
        p->setName("Taxa" + std::to_string(p->getIndex()));
    }
}

/* This method exists to allow a way to quicky get a node via its index. Mainly created to create a traceplot for a given node
   if the index is negative or greater than the number of nodes, this method returns null*/
Node* TreeObject::getNodeWithIndex(int index) {
    if(index < 0 || index > nodes.size()){
        return nullptr;
    }
    for(Node* node :nodes){
        if(node->getIndex() == index){
            return node;
        }
    }
    return nullptr;
}

void TreeObject::flipAllTPs(){
    for(Node* n : nodes){
        n->setNeedsTPUpdate(true);
    }
}

void TreeObject::flipAllCLs(){
    for(Node* n : nodes){
        if(!n->getIsTip()){
            n->setNeedsCLUpdate(true);
        }
    }
}

//For outputting a newick string
void TreeObject::writeNode(Node* p, std::stringstream& strm) const{
    if(p == nullptr)
        return;
    
    if(!p->getIsTip())
        strm << "(";
    else
        strm << p->getName(); 
        // << "[&index=" << p->getIndex() << "]";

    std::set<Node*>& pDesc = p->getNeighbors();
    bool foundFirst = false;
    for(Node* n : pDesc){
        if(n != p->getAncestor()){
            if(foundFirst)
                strm << ",";
            foundFirst = true;
            writeNode(n, strm);
        }
    }

    if(!p->getIsTip())
        strm << ")";
    // if(!p->getIsTip())
    //     strm << ")" << p->getName() <<"[&index=" << p->getIndex() << "]";
    if(p->getAncestor() != nullptr)
        strm << ":" << branchLengths[p->getIndex()];
    else
        strm << ":0.1;";
}


