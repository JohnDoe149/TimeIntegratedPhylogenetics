#include "RateMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>

RateMatrix::RateMatrix(Settings settings) : 
                                   currentQMatrix(4, 4, 0.0), oldQMatrix(4, 4, 0.0), 
                                   currentStationary(4, -1), oldStationary(4, -1), stationaryAlpha(50) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i != j)
                currentQMatrix(i,j) = 1;
        }
    }

    std::vector<double> alpha;
    for(int i = 0; i < 4; i++)
        alpha.push_back(1.0);
    
    Probability::Dirichlet::rv(&rng, alpha, currentStationary);
    oldStationary = currentStationary;
    
    oldQMatrix = currentQMatrix.copy();

    dirty();
}

void RateMatrix::accept() {
    oldQMatrix = currentQMatrix.copy();

    oldStationary = currentStationary;
}

void RateMatrix::reject() {
    currentQMatrix = oldQMatrix.copy();

    currentStationary = oldStationary;
}

double RateMatrix::lnPrior() {
    return 0;
}

double RateMatrix::updateRates(){
    //...
}

double RateMatrix::updateStationary(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    dirty();

    std::vector<double> z(currentStationary.size(), 0.0);
    std::vector<double> alphaForward(currentStationary.size(), 0.0);
    std::vector<double> alphaReverse(currentStationary.size(), 0.0);

    for(int i = 0; i < 4; i++){
        alphaForward[i] = (currentStationary[i] * stationaryAlpha) + 0.1;
    }

    Probability::Dirichlet::rv(&rng, alphaForward, z);

    for(int i = 0; i < 4; i++){
        alphaReverse[i] = (z[i] * stationaryAlpha) + 0.1;
    }

    double hastings = Probability::Dirichlet::lnPdf(alphaReverse, currentStationary) - Probability::Dirichlet::lnPdf(alphaForward, z);

    currentStationary = z;

    return hastings;
}

Matrix<double> RateMatrix::Q() {
    Matrix<double> returnMatrix(currentQMatrix.copy());

    double scaler= 0.0;
    for(int i = 0; i < 4; i++){
        double total = 0.0;
        for(int j = 0; j < 4; j++){
            if(j != i){
                total += returnMatrix(i , j);
            }
        }
        returnMatrix(i, i) = total * -1;
        scaler += returnMatrix(i, i) * currentStationary[i];
    }	

    scaler = -1.0 / scaler;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            returnMatrix(i, j) *= scaler;

    return returnMatrix;
}

void RateMatrix::tune(){}