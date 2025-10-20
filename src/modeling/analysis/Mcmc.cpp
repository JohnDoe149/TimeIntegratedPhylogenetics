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
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;
 
    branchChoice = s.branchWeight;
    topologyChoice = branchChoice + s.topologyWeight;
    stationaryChoice = topologyChoice + s.stationaryWeight; 
    rateChoice = stationaryChoice + s.rateWeight;

    model->regenerateLikelihood();
    model->accept();
}

void Mcmc::burnin(){

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();
    auto gibbsSampler = [this, &rng] (std::function<double()> updater, double currentLnPosterior) {
        double lnProposalRatio = updater();
        model->regenerateLikelihood();

        double modelPrior = model->lnPrior();
        double modelLikelihood = model->lnLikelihood();

        double newLnPosterior = modelLikelihood + modelPrior;
        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + lnPosteriorRatio;
        if(std::log(rng.uniformRv()) < lnR){
            model->accept();
            return newLnPosterior;
        }
        else{
            model->reject();
            return currentLnPosterior;
        }
    };

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        double randomMove = rng.uniformRv() * rateChoice;
        int gibbsUpdates;
        int numNodes = tree->getTree()->getNumNodes();

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
    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();
    auto gibbsSampler = [this, &rng] (std::function<double()> updater, double currentLnPosterior) {
        double lnProposalRatio = updater();
        model->regenerateLikelihood();

        double modelPrior = model->lnPrior();
        double modelLikelihood = model->lnLikelihood();

        double newLnPosterior = modelLikelihood + modelPrior;
        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + lnPosteriorRatio;
        if(std::log(rng.uniformRv()) < lnR){
            model->accept();
            return newLnPosterior;
        }
        else{
            model->reject();
            return currentLnPosterior;
        }
    };

    std::string tabularHeader = model->tabularHeader();
    std::cout << tabularHeader;

    std::ofstream fs;
    fs.open(analysisLog, std::ofstream::out);
    fs << tabularHeader;
    fs.close();

    fs.open(treeLog, std::ofstream::out);
    fs << model->treeHeader();
    fs.close();

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
