#include "core/RandomVariable.hpp"
#include "core/Alignment.hpp"
#include "core/Settings.hpp"
#include "core/Probability.hpp"
#include "ncl/nxscharactersblock.h"
#include "modeling/parameters/trees/TreeObject.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "modeling/parameters/RateMatrix.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/model/TransitionProbability.hpp"
#include "modeling/analysis/Mcmc.hpp"
#include "modeling/parameters/trees/Node.hpp"
#include <algorithm>
#include <chrono>
#include "test.h"

#ifndef TEST
int main(int argc, char* argv[]) {

    Settings settings(argc, argv);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Alignment aln(settings.nexusInput);
    std::cout << "Initializing model..." << std::endl;

    TreeParameter treeParam(&aln, settings.fixedTree, settings.treeLengthLambda);

    RateMatrix rateMatrix(settings);

    Model model(settings, &aln, &treeParam, &rateMatrix);

    Mcmc myMCMC(&model, &treeParam, &rateMatrix, settings);

    std::cout << "Starting MCMC..." << std::endl;

    myMCMC.burnin();
    myMCMC.run();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << treeParam.writeNewick() << "\n" << std::flush;
    std::cout << "Analysis was completed in " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << "[m]" << std::endl;
}
#endif

#ifdef TEST
int main(int argc, char* argv[]) {
    Settings settings(argc, argv);
    RandomVariable& rng = RandomVariable::randomVariableInstance(100);
    Alignment aln(settings.nexusInput);
    TreeParameter treeParam(&aln, settings.fixedTree, settings.treeLengthLambda);
    RateMatrix rateMatrix(settings);
    TreeObject *treeObject = treeParam.getTree();
    TransitionProbability *transProb = new TransitionProbability(treeObject->getNumNodes());

    // // hijack the Q update with a matrix that is guaranteed to be complex
    // Matrix<double> Q(4,4,0.0);
    // Q(0,1) = -1;
    // Q(1,0) = 1;
    // Q(2,2) = 2;
    // Q(2,3) = -3;
    // Q(3,2) = 3;
    // Q(3,3) = 2;
    // Q.print();
    // transProb->updateQ(Q);
    // const std::vector<Node*> poSeq = treeObject->getPostOrderSeq();
    // for(Node* n : poSeq){
    //     int nIndex = n->getIndex();
    //     if(n != treeObject->getRoot()) {
    //         n->setNeedsTPUpdate(true); 
    //         std::vector<double> gammaVec = treeObject->getGammaParams(n);

    //         // inside the setProbs method, we will be able to see the transition matrix p0
    //         transProb->setProbs(0, 0, nIndex, gammaVec[0], gammaVec[1]);
    //     }
    //     n->setNeedsTPUpdate(false);
    //     break;
    // }

    //force updates to stationary and rate matrix
    const std::vector<Node*> poSeq = treeObject->getPostOrderSeq();
    for(int i = 0; i < 30; i ++){
        Matrix<double> Q = rateMatrix.Q();
        transProb->updateQ(Q);
        for(Node* n : poSeq){
            int nIndex = n->getIndex();
            if(n != treeObject->getRoot()) {
                n->setNeedsTPUpdate(true); 
                std::vector<double> gammaVec = treeObject->getGammaParams(n);
                transProb->setProbs(0, 0, nIndex, gammaVec[0], gammaVec[1]);
            }
            n->setNeedsTPUpdate(false);
            break;
        }
    }
}
#endif