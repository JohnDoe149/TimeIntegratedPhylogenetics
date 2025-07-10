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
        RateMatrix(Settings settings);
        RateMatrix();
        Matrix<double> Q();
        void accept();
        void reject();
        void tune();
        double lnPrior();
        double updateRates();
        double updateStationary();
        std::vector<double> getStationary() {return currentStationary;}
        std::vector<double> getRate() {return nucleotideTransitionRates;}
        std::vector<std::string> transNameOrder(){return {"AC", "AG", "AT", "CG", "CT", "GT"};} 
        int rateAcceptCount;
        int rateCount;
        int stationaryAcceptCount;
        int stationaryCount;

        // if 1, rate was selected, if 2, stationary was selected
        int rateOrStationary = 0;
    private:
        std::vector<double> nucleotideTransitionRates;
        std::vector<double> oldNucleotideTransitionRates;

        double stationaryAlpha;
        double rateStepsize;
        double stationaryStepsize;

        std::vector<double> currentStationary;
        std::vector<double> oldStationary;

        double ratePrior;
        double oldRatePrior;
        double stationaryPrior;
        double oldStationaryPrior;
};

#endif