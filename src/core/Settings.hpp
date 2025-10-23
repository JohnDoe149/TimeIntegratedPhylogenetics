#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#include <string>
#include <vector>

struct Settings {
                    Settings(void) = delete;
                    Settings(int argc,  char* argv[]); // initialize settings from command line arguments or set manually in class constructor
    std::string     nexusInput; // path to MSA input
    std::string     treeOutput; // path to tree log output
    std::string     mcmcOutput; // path to MCMC trace output

    std::string     fixedTree; // NEWICK string for fixed tree analyses, empty string if no fixed tree is used
    int             numIterations; // number of MCMC iterations, does not include burn-in
    int             printFrequency; // frequency to print to terminal
    int             sampleFrequency; // frequency to log to file
    int             burnInIterations; // number of burn-in iterations
    int             tuneFrequency; // frequency to tune moves during burn-in
    double          treeLengthLambda; // lambda parameter for tree length prior
    double          branchWeight; // weight for branch length moves during MCMC (gamma branch params)
    double          topologyWeight; // weight for topology moves during MCMC (NNI)
    double          stationaryWeight; // weight for stationary distribution moves during MCMC
    double          rateWeight; // weight fortransition rate moves during MCMC

    void            usage(); // prints out to terminal all possible settings
    void            print(); // prints out to terminal what all settings are set to
};

#endif
