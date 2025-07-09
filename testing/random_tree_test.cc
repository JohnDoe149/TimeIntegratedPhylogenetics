#include <gtest/gtest.h>
#include "modeling/model/RandomTree.hpp"
#include "modeling/parameters/trees/Node.hpp"

// verify that the preorder traversal works
// given that right and left arent fixed, as long as pre order and post order are roughly in reverse order preorder is right
TEST(RandomTree, preorderTraversal){
    RandomTree randomTree(4, 10);
    std::vector<Node*> preorder = randomTree.getPreorderTrav();
    std::vector<Node*> postorder = randomTree.getTree()->getPostOrderSeq();
    // for(int i = 0; i < preorder.size(); i++){
    //     std::cout << preorder[i] << "pre has an index of " << (*preorder[i]).getIndex()<< "\n";
    // }
    // for(int i = 0; i < preorder.size(); i++){
    //     std::cout << postorder[i] << "post has an index of " << (*postorder[i]).getIndex()<< "\n";
    // }
}

// verify that currentSequences look good when genNewData is called
TEST(RandomTree, genNewData){
    RandomTree randomTree(4, 10);
    randomTree.genNewData();
    std::vector<std::vector<int>> allSeq = randomTree.getAllNodeSequences();
    for(int i = 0; i < allSeq.size(); i++){
        for(int j = 0; j < allSeq[i].size(); i++){
            std::cout << allSeq[i][j];
        }
        std::cout << "\n";
    }
}