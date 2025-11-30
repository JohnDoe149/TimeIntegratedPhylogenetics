#include "Mcmc.hpp"
#include "core/RandomVariable.hpp"
#include "modeling/model/Model.hpp"
#include "modeling/parameters/Parameter.hpp"
#include "modeling/parameters/RateMatrix.hpp"
#include "modeling/parameters/trees/TreeParameter.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include "test.h"

Mcmc::Mcmc(Model* m, TreeParameter* t, RateMatrix* rm, Settings& s) : 
    model(m), rateMatrix(rm), tree(t) { 

    // set MCMC parameters from settings
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
        
    // set move weights from settings
    branchChoice = s.branchWeight;
    topologyChoice = branchChoice + s.topologyWeight;
    stationaryChoice = topologyChoice + s.stationaryWeight; 
    rateChoice = stationaryChoice + s.rateWeight;

    // initialize MCMC by calculating initial likelihood
    model->regenerateLikelihood();
    model->accept();
}

void Mcmc::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior(); // get the last accepted posterior

    // define a gibbs sampler lambda function to update parameters
    auto gibbsSampler = [this, &rng] (std::function<double()> updater, double currentLnPosterior) {
        double lnProposalRatio = updater(); // call updater to propose a new model state, returns the hasting ratio
        model->regenerateLikelihood(); // recalculate likelihood after proposal

        double modelPrior = model->lnPrior(); // get the prior after proposal
        double modelLikelihood = model->lnLikelihood(); // get the new likelihood
        double newLnPosterior = modelLikelihood + modelPrior; // calculate new posterior

        double lnPosteriorRatio = newLnPosterior - currentLnPosterior; // calculate posterior ratio
        double lnR = lnProposalRatio + lnPosteriorRatio; // apply the hastings ratio
        if(std::log(rng.uniformRv()) < lnR){ // acceptance check based on a "dice roll"
            model->accept();
            return newLnPosterior;
        }
        else{
            model->reject();
            return currentLnPosterior;
        }
    };

    // begin burn-in iterations
    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }


        double randomMove = rng.uniformRv() * rateChoice; // randomly pick a move based on weights
        int gibbsUpdates;
        int numNodes = tree->getTree()->getNumNodes();

        // Employ a blocked gibbs sampler, for topology and branch updates
        if (randomMove < topologyChoice){
            std::function<double()> topologyUpdater = [this]() { return tree->updateTreeMove(); };
            std::function<double()> branchUpdater = [this]() { return tree->updateTreeGamma(); };
            gibbsUpdates = numNodes / 2;
            for (int i = 0; i < gibbsUpdates; i++){
                currentLnPosterior = gibbsSampler(topologyUpdater, currentLnPosterior);
                currentLnPosterior = gibbsSampler(branchUpdater, currentLnPosterior);
                currentLnPosterior = gibbsSampler(branchUpdater, currentLnPosterior);
            }
        }
        else if(randomMove < stationaryChoice){
            std::function<double()> stationaryUpdater = [this]() { return rateMatrix->updateStationary(); };
            gibbsUpdates = 2;
            for (int i = 0; i < gibbsUpdates; i++){
                currentLnPosterior = gibbsSampler(stationaryUpdater, currentLnPosterior);
            }
        }
        else {
            std::function<double()> rateUpdater = [this]() { return rateMatrix->updateRates(); };
            gibbsUpdates = 10;
            for (int i = 0; i < gibbsUpdates; i++){
                currentLnPosterior = gibbsSampler(rateUpdater, currentLnPosterior);
            }
        }

    }
}

void Mcmc::run(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior(); // get the last accepted posterior

    // define a gibbs sampler lambda function to update parameters
    auto gibbsSampler = [this, &rng] (std::function<double()> updater, double currentLnPosterior) {
        double lnProposalRatio = updater(); // call updater to propose a new model state, returns the hasting ratio
        model->regenerateLikelihood(); // recalculate likelihood after proposal

        double modelPrior = model->lnPrior(); // get the prior after proposal
        double modelLikelihood = model->lnLikelihood(); // get the new likelihood

        double newLnPosterior = modelLikelihood + modelPrior; // calculate new posterior
        double lnPosteriorRatio = newLnPosterior - currentLnPosterior; // calculate posterior ratio
        double lnR = lnProposalRatio + lnPosteriorRatio; // apply the hastings ratio
        if(std::log(rng.uniformRv()) < lnR){ // acceptance check based on a "dice roll"
            model->accept();
            return newLnPosterior;
        }
        else{
            model->reject();
            return currentLnPosterior;
        }
    };

    // sort logging output
    std::string tabularHeader = model->tabularHeader();
    std::cout << tabularHeader;

    std::ofstream fs;
    fs.open(analysisLog, std::ofstream::out);
    fs << tabularHeader;
    fs.close();

    fs.open(treeLog, std::ofstream::out);
    fs << model->treeHeader();
    fs.close();

    // begin sampling iterations
    for(int n = 1; n <= numIter; n++){
        if(n % printFreq == 0){
            std::cout << model->tabularOut(n);
        }
        if(n % sampleFreq == 0){
            fs.open(analysisLog, std::ofstream::app);
            fs << model->tabularOut(n);
            fs.close();
            fs.clear();

            fs.open(treeLog, std::ofstream::app);
            fs << model->treeOut(n);
            fs.close();
            fs.clear();
        }

        double randomMove = rng.uniformRv() * rateChoice;
        int gibbsUpdates;
        int numNodes = tree->getTree()->getNumNodes();

        // Employ a blocked gibbs sampler, for topology and branch updates
        if (randomMove < topologyChoice){
            std::function<double()> topologyUpdater = [this]() { return tree->updateTreeMove(); };
            std::function<double()> branchUpdater = [this]() { return tree->updateTreeGamma(); };
            gibbsUpdates = numNodes / 2;
            for (int i = 0; i < gibbsUpdates; i++){
                currentLnPosterior = gibbsSampler(topologyUpdater, currentLnPosterior);
                currentLnPosterior = gibbsSampler(branchUpdater, currentLnPosterior);
                currentLnPosterior = gibbsSampler(branchUpdater, currentLnPosterior);
            }
        }
        else if(randomMove < stationaryChoice){
            std::function<double()> stationaryUpdater = [this]() { return rateMatrix->updateStationary(); };
            gibbsUpdates = 2;
            for (int i = 0; i < gibbsUpdates; i++){
                currentLnPosterior = gibbsSampler(stationaryUpdater, currentLnPosterior);
            }
        }
        else {
            std::function<double()> rateUpdater = [this]() { return rateMatrix->updateRates(); };
            gibbsUpdates = 10;
            for (int i = 0; i < gibbsUpdates; i++){
                currentLnPosterior = gibbsSampler(rateUpdater, currentLnPosterior);
            }
        }
    }
}
