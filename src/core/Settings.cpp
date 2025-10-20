#include "Settings.hpp"
#include "Msg.hpp"
#include <iostream>
#include <string>
#include <vector>

Settings::Settings(int argc,  char* argv[]) : nexusInput(""), treeOutput(""), mcmcOutput(""),
                                              numIterations(100000), printFrequency(10), sampleFrequency(10),
                                              burnInIterations(10000), tuneFrequency(100), treeLengthLambda(1.0),
                                              stationaryWeight(6.0), branchWeight(15.0), topologyWeight(15.0), rateWeight(8.0), fixedTree("") {

    std::vector<std::string> settings;
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i]; 
        settings.push_back(arg);
    }

    settings.push_back("-nexus");
    settings.push_back("/workspaces/JohnPhylo/validation/TIP_testing_data/taxa50diam0.250000siteCount800/0tipSequence.fasta");
    // settings.push_back("/workspaces/JohnPhylo/res/replicase.nex");
    settings.push_back("-treeOut");
    settings.push_back("/workspaces/JohnPhylo/validation/mixTree.tree");
    // settings.push_back("/workspaces/JohnPhylo/res/trees.trees");
    settings.push_back("-mcmcOut");
    settings.push_back("/workspaces/JohnPhylo/validation/mixTree.log");
    // settings.push_back("/workspaces/JohnPhylo/res/analysis.log");

    if (settings.size() == 0) {
        usage();
        Msg::error("Expected command line arguments");
    }

    std::string currentArg = "";
    for (int i=0; i<settings.size(); i++) {
        if (currentArg == "")
            currentArg = settings[i];
        else {
            if (currentArg == "-nexus")
                nexusInput = settings[i];
            else if (currentArg == "-treeOut")
                treeOutput = settings[i];
            else if (currentArg == "-mcmcOut")
                mcmcOutput = settings[i];
            else if (currentArg == "-numIter")
                numIterations = stoi(settings[i]);
            else if (currentArg == "-printFreq")
                printFrequency = stoi(settings[i]);
            else if (currentArg == "-sampleFreq")
                sampleFrequency = stoi(settings[i]);
            else if (currentArg == "-burnInIter")
                burnInIterations = stoi(settings[i]);
            else if (currentArg == "-tuneFreq")
                tuneFrequency = stoi(settings[i]);
            else if (currentArg == "-treeLambda")
                treeLengthLambda = stod(settings[i]);
            else if (currentArg == "-stationaryWeight")
                stationaryWeight = stod(settings[i]);
            else if (currentArg == "-rateWeight")
                rateWeight = stod(settings[i]);
            else if (currentArg == "-branchWeight")
                branchWeight = stod(settings[i]);
            else if (currentArg == "-topologyWeight")
                topologyWeight = stod(settings[i]);
            else if (currentArg == "-fixedTree")
                fixedTree = settings[i];
            else{
                Msg::error("Could not interpret argument " + currentArg);
                usage();
            }
            currentArg = "";
        }
    }

    if(nexusInput == "" || treeOutput == "" || mcmcOutput == ""){
        usage();
        Msg::error("For non-simulation analyses, nexusInput, treeOutput, and mcmcOutput are required arguments.");
    }

    print();
}


void Settings::print(){
    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : " << nexusInput << std::endl;
    std::cout << "   * -treeOut           : " << treeOutput << std::endl;
    std::cout << "   * -mcmcOut           : " << mcmcOutput << std::endl;
    std::cout << "   * -fixedTree         : " << fixedTree << std::endl;
    std::cout << std::endl;
    
    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : " << treeLengthLambda << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : " << numIterations << std::endl;
    std::cout << "   * -printFreq         : " << printFrequency << std::endl;
    std::cout << "   * -sampleFreq        : " << sampleFrequency << std::endl;
    std::cout << "   * -burnInIter        : " << burnInIterations << std::endl;
    std::cout << "   * -tuneFreq          : " << tuneFrequency << std::endl;
    std::cout << "   * -topologyWeight    : " << topologyWeight << std::endl;
    std::cout << "   * -branchWeight      : " << branchWeight << std::endl;
    std::cout << "   * -stationaryWeight  : " << stationaryWeight << std::endl;
    std::cout << "   * -rateWeight        : " << rateWeight << std::endl;
    std::cout << std::endl;
}

void Settings::usage(void) {

    std::cout << "Inference Input/Output:" << std::endl;
    std::cout << "   * -nexus             : Input nexus file containing the nculeotide alignment." << std::endl;
    std::cout << "   * -treeOut           : The output file name for the tree trace." << std::endl;
    std::cout << "   * -mcmcOut           : The output file name for the bulk of the MCMC trace." << std::endl;
    std::cout << "   * -fixedTree         : The NEWICK string corresponding to the fixed tree you wish to analyze." << std::endl;
    std::cout << std::endl;

    std::cout << "Model Parameters:" << std::endl;
    std::cout << "   * -treeLambda        : Lambda parameter for the tree length exponential prior." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Sampling Options:" << std::endl;
    std::cout << "   * -numIter           : The number of iterations for the MCMC." << std::endl;
    std::cout << "   * -printFreq         : How often to output the MCMC state to the screen." << std::endl;
    std::cout << "   * -sampleFreq        : How often to ouput the MCMC state to log files." << std::endl;
    std::cout << "   * -burnInIter        : The number of iterations for the burn-in." << std::endl;
    std::cout << "   * -tuneFreq          : How often to tune the MCMC moves during the burn-in." << std::endl;
    std::cout << "   * -treeWeight        : How often to propose a move on the tree." << std::endl;
    std::cout << "   * -stationaryWeight  : How often to propose a move on the stationary distribution." << std::endl;
    std::cout << "   * -rateWeight        : How often to propose a move on the rate parameters." << std::endl;
    std::cout << std::endl;
}