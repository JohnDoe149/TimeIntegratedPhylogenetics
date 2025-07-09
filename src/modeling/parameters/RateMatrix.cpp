#include "RateMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include "test.h"

RateMatrix::RateMatrix(Settings settings) : 
                                   currentQMatrix(4, 4, 0.0), oldQMatrix(4, 4, 0.0), 
                                   currentStationary(4, -1), oldStationary(4, -1), stationaryAlpha(50), rateStepsize(0.3), 
                                   rateAcceptCount(0), rateCount(0), stationaryAcceptCount(0), stationaryCount(0){

    #ifdef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(100);
    #endif 
    #ifndef TEST
    RandomVariable& rng = RandomVariable::randomVariableInstance(); 
    #endif 
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i != j)
                currentQMatrix(i,j) = 1;
        }
    }

    // fill a vector alpha with all ones and generate a random draw of stationaries that sum to 1
    std::vector<double> alpha;
    for(int i = 0; i < 4; i++)
        alpha.push_back(2.0);
    Probability::Dirichlet::rv(&rng, alpha, currentStationary);
    oldStationary = currentStationary;
    oldQMatrix = currentQMatrix.copy();
    dirty();
}

RateMatrix::RateMatrix() : 
                            currentQMatrix(4, 4, 0.0), oldQMatrix(4, 4, 0.0), 
                            currentStationary(4, -1), oldStationary(4, -1), stationaryAlpha(50), rateStepsize(0.3), 
                            rateAcceptCount(0), rateCount(0), stationaryAcceptCount(0), stationaryCount(0){
    RandomVariable& rng = RandomVariable::randomVariableInstance(); 
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i != j)
                currentQMatrix(i,j) = 1;
        }
    }

    // fill a vector alpha with all ones and generate a random draw of stationaries that sum to 1
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
    if(rateOrStationary == 2){
        rateAcceptCount++;
    }
    else if (rateOrStationary == 1){
        stationaryAcceptCount++;
    }
}

void RateMatrix::reject() {
    currentQMatrix = oldQMatrix.copy();
    currentStationary = oldStationary;
    rateOrStationary = 0;
}

double RateMatrix::lnPrior() {
    return 0;
}

// make an update to the rates
double RateMatrix::updateRates(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    rateCount++;
    dirty();
    rateOrStationary = 2; // set this to 2 to indicate that the choice has been made to update rate

    std::vector<double> alpha;
    // this vector will store all of our rates, takes the upper right values left to right row major order
    for(int i = 0; i<=3; i++){
        for(int j = i+1; j<=3; j++){
            alpha.push_back(currentQMatrix(i,j));
        }
    }

    double log_hastings = 0;
    // make some proposals for the new rates
    // we need to propose individual rates as the rates are not correlated
    for(int i = 0; i<alpha.size(); i++){
        // scaling factor (stepsize)
        double scale = std::exp(rateStepsize * (rng.uniformRv() - 0.5));
        alpha[i] = alpha[i] * scale;
        log_hastings += std::log(scale);
    }
    
    // now slide the new vectors back into the new matrix
    int indexer = 0;
    for(int i = 0; i<=3; i++){
        for(int j = i+1; j<=3; j++){
           double newRate = alpha[indexer];
           currentQMatrix(i, j) = newRate;
           currentQMatrix(j, i) = newRate;
           indexer++;
        }
    }
    currentQMatrix = Q();
    return log_hastings;
}

double RateMatrix::updateStationary(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    dirty();
    stationaryCount++;
    rateOrStationary = 1; // set this to 1 to indicate that stationaries are being updated

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

            // if not diagonal, sum up total so we can make the diagonal the additive inverse
            // returnMatrix.print();
            // std::cout << "break \n";
            if(j != i){
                returnMatrix(i,j) = currentStationary[j] * returnMatrix(i,j);
                total += returnMatrix(i , j);
            } 
            // returnMatrix.print();
            // std::cout << "break \n";
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

void RateMatrix::tune(){ 
    double transitionRateRate = (double) rateAcceptCount/rateCount;
    if(transitionRateRate > 0.33){
        rateStepsize *= (1.0 + ((transitionRateRate-0.33)/0.67));
    }
    else {
        rateStepsize /= (2.0 - transitionRateRate/0.33);
    }
    rateAcceptCount = 0;
    rateCount = 0;

    // we need to modify the stationaryAlpha
    double stationaryRate = (double) stationaryAcceptCount/stationaryCount;
    // if our acceptance rate is too high, that means we are exploring too much, so dial it back a little bit
    if(stationaryRate > 0.33){
        stationaryAlpha /= (1.0 + ((stationaryRate-0.33)/0.67));
    }
    else {
        stationaryAlpha *= (2.0 + ((stationaryRate-0.33)));
    }
    stationaryAcceptCount = 0;
    stationaryCount = 0;
}