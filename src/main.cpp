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
#include "modeling/model/RandomTree.hpp"

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

    RandomTree randomTree(6, 100);
    randomTree.genNewData();
    std::vector<std::vector<int>> allSeq = randomTree.getAllNodeSequences();
    for(int i = 0; i < allSeq.size(); i++){
        for(int j = 0; j < allSeq[i].size(); j++){
            std::cout << allSeq[i][j];
        }
        std::cout << "\n";
    }
}
#endif