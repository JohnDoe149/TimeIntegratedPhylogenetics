#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include <string>

class Node{
    public:
                            Node();
        void                addNeighbor(Node* n) {neighbors.insert(n);} // adds n to its neighbor set
        static Node*        chooseNodeFromSet(std::set<Node*>& s); // given a set of nodes, return a random node from the set
        Node*               getAncestor() {return ancestor;} // return the ancestor of the node
        void                setAncestor(Node* a) {ancestor = a;} // set the ancestor of the node
        int                 getIndex() const {return index;} // return the unique index of the node
        bool                getIsTip() const {return isTip;} // return whether the node is a tip/leaf
        std::string         getName() const {return name;} // return the name of the node
        bool                getNeedsTPUpdate(){return needsTPUpdate;} // return whether node needs transition probability update
        std::set<Node*>&    getNeighbors() {return neighbors;} // return the set of neighboring nodes
        int                 getOffset() const {return offset;} // return the offset
        void                removeNeighbor(Node* n) {neighbors.erase(n);} // removes n from its neighbor set, does nothing if n is absent
        void                removeAllNeighbors() {neighbors.clear();} // removes all neighbors its neighbor set
        void                setIndex(int i) {index = i;} // set the unique index for the node
        void                setIsTip(bool t) {isTip = t;} // set whether node is a tip/leaf
        void                setName(std::string s) {name = s;} // set the name of the node
        void                setNeedsTPUpdate(bool nU){needsTPUpdate = nU;} // set whether node needs transition probability update
        void                setOffset(int o){offset = o;} // update node offset for new nodes added
    private:
        Node*               ancestor; // all nodes have an ancestor except the root
        int                 index; // unique index for the node
        bool                isTip; // whether the node is a tip/leaf
        std::string         name; // name of the node, all tips/leaves have names from the alignment
        bool                needsTPUpdate; // whether the node needs its transition probabilities updated
        std::set<Node*>     neighbors; // set of neighboring nodes, includes ancestors and children
        int                 offset; // used during tree construction for identification
};

#endif