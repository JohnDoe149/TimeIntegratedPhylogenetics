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
#include <filesystem>
#include <iomanip>

#include "test.h"
#include "modeling/model/RandomTree.hpp"

#ifdef BASE
int main(int argc, char* argv[]) {

    // command line arguments can be set ahead of time in settings.cpp
    Settings settings(argc, argv);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Alignment aln(settings.nexusInput, 1);
    std::cout << "Initializing model..." << std::endl;

    TreeParameter treeParam(&aln);

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

#ifdef GENDATA
int main(int argc, char* argv[]) {
    std::vector<double> diameters = {0.25, 0.5, 1.0, 2.0, 3.0, 4.0};
    for (double treeDiameter: diameters){
        int taxa = 50;
        int siteCount = 800;
        double branchLengthLambda = 10.0;
        std::filesystem::path projectRoot = std::filesystem::current_path();
        projectRoot = projectRoot.parent_path();
        std::filesystem::path outputPath = projectRoot / "validation" / "TIP_testing_data" / ("taxa" + std::to_string(taxa) + "diam" + std::to_string(treeDiameter) + "siteCount" + std::to_string(siteCount));

        // make a file path for tree file
        for(int dataCount = 0; dataCount < 5; dataCount++){
            std::ofstream treeOutFile((outputPath / (std::to_string(dataCount) + "newick.tree")));

            // we must set the file path
            RandomTree randomTree(taxa, siteCount, branchLengthLambda, outputPath);
            randomTree.rescaleTree(treeDiameter);
            randomTree.genNewData();
            randomTree.setOutputPath(outputPath);
            TreeObject* randomtree = randomTree.getTree();

            treeOutFile << randomtree->getNewick() << std::endl;

            std::string sequenceName = std::to_string(dataCount) + "completeSequence.fasta";
            std::string tipsName = std::to_string(dataCount) + "tipSequence.fasta";


        }
    }
}
#endif