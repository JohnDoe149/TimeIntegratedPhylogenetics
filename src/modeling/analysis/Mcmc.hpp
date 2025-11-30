#ifndef MCMC_HPP
#define MCMC_HPP
#include <vector>
#include <string>

class Model;
class Parameter;
class TreeParameter;
class RateMatrix;
class Settings;

// Mcmc class serves as a controller for running the model. It handles class calls for burn-in, sampling, tuning, and logging.
class Mcmc{
    public:
        Mcmc(void)=delete;
        Mcmc(Model* m, TreeParameter* t, RateMatrix* rm, Settings& s); // initialize Mcmc with model, treeParameter, rateMatrix as parameters and settings for settings like numIterations, burnInIterations, ETC
        void burnin(); // performs burn-in by running numBurnIn iterations, tuning moves every tuneFreq iterations
        void run(); // performs sampling by running numIter iterations, printing to screen every printFreq iterations and logging every sampleFreq iterations
    private:
        int numIter; // number of iterations for run (does not include burn-in)
        int numBurnIn; // number of iterations for burn-in
        int printFreq; // frequency to print to screen
        int tuneFreq; // frequency to tune moves during burn-in
        int sampleFreq; // frequency to log samples during run to tree and trace files

        TreeParameter* tree; // treeParameter used in the model
        RateMatrix* rateMatrix; // rateMatrix used in the model

        std::string analysisLog; // string to hold the analysis log that will be written to file, check settings/terminal input for file location
        std::string treeLog; // string to hold the tree log that will be written to file, check settings/terminal input for file location

        double branchChoice; // weight for choosing branch move during MCMC
        double topologyChoice; // weight for choosing topology move during MCMC
        double stationaryChoice; // weight for choosing stationary move during MCMC
        double rateChoice; // weights for choosing stationary and rate moves during MCMC

        Model* model; // used to do likelihood calculations and get lnPrior
};

#endif