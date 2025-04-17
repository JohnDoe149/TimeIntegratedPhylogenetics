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
    private:
        Matrix<double> currentQMatrix;
        Matrix<double> oldQMatrix;

        double stationaryAlpha;

        std::vector<double> currentStationary;
        std::vector<double> oldStationary;
};

#endif