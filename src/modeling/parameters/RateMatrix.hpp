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
        Matrix<double> Q();
        void accept();
        void reject();
        void tune();
        double lnPrior();
        double updateRates();
        double updateStationary();
        std::vector<double> getStationary() {return currentStationary;}
        Matrix<double> getRate() {return currentQMatrix;}
        std::vector<std::string> transNameOrder(){return {"AC", "AG", "AT", "CG", "CT", "GT"};}
        int rateAcceptCount;
        int rateCount;
        int stationaryAcceptCount;
        int stationaryCount;
        // if 1, rate was selected, if 2, stationary was selected
        int rateOrStationary = 0;
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        double stationaryAlpha;
        double rateStepsize;
        double stationaryStepsize;

        std::vector<double> currentStationary;
        std::vector<double> oldStationary;
};

#endif