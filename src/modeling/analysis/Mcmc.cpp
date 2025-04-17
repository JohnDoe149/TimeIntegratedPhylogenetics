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

Mcmc::Mcmc(Model* m, TreeParameter* t, RateMatrix* rm, Settings& s) : 
    model(m), rateMatrix(rm), tree(t) { 
    numIter = s.numIterations;
    numBurnIn = s.burnInIterations;
    printFreq = s.printFrequency;
    tuneFreq = s.tuneFrequency;
    sampleFreq = s.sampleFrequency;
    analysisLog = s.mcmcOutput;
    treeLog = s.treeOutput;

    treeChoice = s.treeWeight;
    stationaryChoice = treeChoice + s.stationaryWeight;
    rateChoice = stationaryChoice + s.rateWeight;

    model->regenerateLikelihood();
    model->accept();
}

void Mcmc::burnin(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

    for(int n = 1; n <= numBurnIn; n++){
        if(n % printFreq == 0){
            std::cout << "Burn-in Iteration " << n << ": " << currentLnPosterior << std::endl;
        }
        if(n % tuneFreq == 0){
            model->tuneMoves();
        }

        double randomMove = rng.uniformRv() * rateChoice;

        std::function<double()> updater;

        if(randomMove < treeChoice){
            updater = [this]() { return tree->update(); };
        }
        else if(randomMove < stationaryChoice){
            updater = [this]() { return rateMatrix->updateStationary(); };
        }
        else {
            updater = [this]() { return rateMatrix->updateRates(); };
        }

        double lnProposalRatio = updater();
        model->regenerateLikelihood();

        double newLnPosterior = model->lnLikelihood() + model->lnPrior();

        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + lnPosteriorRatio;

        if(std::log(rng.uniformRv()) < lnR){
            model->accept();
            currentLnPosterior = newLnPosterior;
        }
        else{
            model->reject();
        }
    }
}

void Mcmc::run(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double currentLnPosterior = model->lnLikelihood() + model->lnPrior();

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

        std::function<double()> updater;

        if(randomMove < treeChoice){
            updater = [this]() { return tree->update(); };
        }
        else if(randomMove < stationaryChoice){
            updater = [this]() { return rateMatrix->updateStationary(); };
        }
        else {
            updater = [this]() { return rateMatrix->updateRates(); };
        }

        double lnProposalRatio = updater();
        model->regenerateLikelihood();

        double newLnPosterior = model->lnLikelihood() + model->lnPrior();

        double lnPosteriorRatio = newLnPosterior - currentLnPosterior;
        double lnR = lnProposalRatio + lnPosteriorRatio;

        if(std::log(rng.uniformRv()) < lnR){
            model->accept();
            currentLnPosterior = newLnPosterior;
        }
        else{
            model->reject();
        }
    }
}
