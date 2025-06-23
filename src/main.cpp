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
    Matrix<double> Q = rateMatrix.Q();
    TreeObject *treeObject = treeParam.getTree();
    TransitionProbability *transProb = new TransitionProbability(treeObject->getNumNodes());
    transProb->updateQ(Q);

    const std::vector<Node*> poSeq = treeObject->getPostOrderSeq();
    for(Node* n : poSeq){
        int nIndex = n->getIndex();
        if(n != treeObject->getRoot()) {
            n->setNeedsTPUpdate(true);
            std::vector<double> gammaVec = treeObject->getGammaParams(n);

            // inside the setProbs method, we will be able to see the transition matrix p0
            transProb->setProbs(0, 0, nIndex, gammaVec[0], gammaVec[1]);
            return 1;
        }
        n->setNeedsTPUpdate(false);
    }
}
#endif