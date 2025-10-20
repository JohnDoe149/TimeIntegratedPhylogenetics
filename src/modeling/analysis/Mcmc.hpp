#ifndef MCMC_HPP
#define MCMC_HPP
#include <vector>
#include <string>

class Model;
class Parameter;
class TreeParameter;
class RateMatrix;
class Settings;

class Mcmc{
    public:
        Mcmc(void)=delete;
        Mcmc(Model* m, TreeParameter* t, RateMatrix* rm, Settings& s);
        void burnin();
        void run();
    private:
        int numIter;
        int numBurnIn;
        int printFreq;
        int tuneFreq;
        int sampleFreq;

        TreeParameter* tree;
        RateMatrix* rateMatrix;

        std::string analysisLog;
        std::string treeLog;

        double branchChoice;
        double topologyChoice;
        double stationaryChoice;
        double rateChoice;

        Model* model;
};

#endif