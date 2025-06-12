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

    // I just want to test gamma map and make sure it works
    TreeObject *treeObject = treeParam.getTree();
    std::vector<Node*> nodes = treeObject->getPostOrderSeq();
    std::vector<double> gammaParams;
    gammaParams.reserve(sizeof(double) * 2 );

    for(Node* n: nodes){
        treeObject->setGammaDist(n, 100*rng.uniformRv(), 100*rng.uniformRv());
        gammaParams = treeObject->getGammaParams(n);
    }
    std::vector<std::vector<double>> allGammas = treeObject->getGammas();

}
#endif