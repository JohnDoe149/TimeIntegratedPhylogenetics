#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#include <string>
#include <vector>

struct Settings {
                    Settings(void) = delete;
                    Settings(int argc,  char* argv[]);
    std::string     nexusInput;
    std::string     treeOutput;
    std::string     mcmcOutput;

    std::string     fixedTree;
    int             numIterations;
    int             printFrequency;
    int             sampleFrequency;
    int             burnInIterations;
    int             tuneFrequency;
    double          treeLengthLambda;
    double          branchWeight;
    double          topologyWeight;
    double          stationaryWeight;
    double          rateWeight;

    void            usage();
    void            print();
};

#endif
