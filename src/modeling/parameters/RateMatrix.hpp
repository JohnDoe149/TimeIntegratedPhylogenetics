#ifndef RATE_MATRIX_HPP
#define RATE_MATRIX_HPP
#include "core/Matrix.hpp"
#include "core/Msg.hpp"
#include "modeling/parameters/Parameter.hpp"
#include <set>
#include <vector>

class Settings;

class RateMatrix : public Parameter {
    public:
        RateMatrix(Settings settings); // initialize a rateMatrix based on param settings
        RateMatrix(); // initialize a base rateMatrix
        Matrix<double> Q(); // generates the Q matrix using current t-rates and stationaries, also normalizes
        void accept(); // call to accept either proposed stationary or rate update
        void reject(); // call to reject either proposed stationary or rate update
        void tune(); // call during burn-in to adjust step-size for proposals. stationary and rate have seperate step-sizes
        double lnPrior(); // returns the prior of current stationary + rate
        double updateRates(); // pick a single transition rate within the matrix and propose a rescaling update.
        double updateStationary(); // pick a single stationary and propose a Dirichlet update
        std::vector<double> getStationary() {return currentStationary;} // return vector of stationary frequences (4 frequences for nucleotides)
        std::vector<double> getRate() {return nucleotideTransitionRates;}   // return vector of transition rates (6 rates for GTR)
        std::vector<std::string> transNameOrder(){return {"AC", "AG", "AT", "CG", "CT", "GT"};}  // return transition names for log printout
        int rateAcceptCount; // number of rate updates accepted
        int rateCount; // number of rate updates proposed
        int stationaryAcceptCount; // number of stationary updates accepted
        int stationaryCount; // number of stationary updates proposed
        int rateOrStationary = 0;  // if 1, rate was selected, if 2, stationary was selected, used to track acceptance for tuning during burn in
    private:
        std::vector<double> nucleotideTransitionRates; // vector of 6 transition rates for GTR representing the proposed transition rates
        std::vector<double> oldNucleotideTransitionRates; // vector of 6 transition rates for GTR representing the previous accepted transition rates 

        double stationaryAlpha; // controls step-size for stationary
        double rateStepsize;  // controls step-size for rate

        std::vector<double> currentStationary; // vector of 4 stationaries representing the proposed stationary frequencies
        std::vector<double> oldStationary; // vector of 4 stationaries representing the previous accepted stationaries

        double ratePrior; // prior of the proposed rates
        double oldRatePrior; // prior of the previous accepted rates
        double stationaryPrior; // prior of the proposed stationaries
        double oldStationaryPrior; // prior of the previous accepted stationaries
};

#endif