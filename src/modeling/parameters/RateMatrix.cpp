#include "RateMatrix.hpp"
#include "core/Matrix.hpp"
#include "core/RandomVariable.hpp"
#include "core/Probability.hpp"
#include "core/Settings.hpp"
#include <cmath>
#include <algorithm>
#include "test.h"

/* This class handles the proposal, accept, reject and general maintenance of rateMatrix. A rateMatrix in this model
   is a GTR-model with 6 nucleotide transition rates (A->C, A->G, etc) and 4 stationary distribution values (A=0.26, C = 0.24, ETC)
   NOTE: The actual rateMatrix is not stored as an instance variable but is generated when rateMatrix::Q()is called using the 
   transitionrates and stationary distribution.
*/
RateMatrix::RateMatrix(Settings settings) : 
                                   currentStationary(4, -1), oldStationary(4, -1), stationaryAlpha(50), rateStepsize(0.3), 
                                   rateAcceptCount(0), rateCount(0), stationaryAcceptCount(0), stationaryCount(0),
                                   nucleotideTransitionRates(6.0, 0.0), oldNucleotideTransitionRates(6.0, 0.0){

    // fill a vector alpha with all ones and generate a random draw of stationaries that sum to 1
    RandomVariable& rng = RandomVariable::randomVariableInstance(); 
    std::vector<double> alpha(4, 2.0);
    Probability::Dirichlet::rv(&rng, alpha, currentStationary);
    oldStationary = currentStationary;

    // set the prior for stationary
    stationaryPrior = Probability::Dirichlet::lnPdf(alpha, currentStationary);

    // now set the nucleotideTransitionRates and the ratePrior
    ratePrior = 0;
    for(int nucIter = 0; nucIter < nucleotideTransitionRates.size(); nucIter++){
        nucleotideTransitionRates[nucIter] = Probability::Gamma::rv(&rng, 2, 5);
        ratePrior += Probability::Gamma::lnPdf(2, 5,nucleotideTransitionRates[nucIter]);
    }

    oldStationaryPrior = stationaryPrior;
    oldRatePrior = ratePrior;
    dirty();
}

// a constructor that does not take in settings as a parameter, meant to be used for testing. 
RateMatrix::RateMatrix() : 
                    currentStationary(4, -1), oldStationary(4, -1), stationaryAlpha(50), rateStepsize(0.3), 
                    rateAcceptCount(0), rateCount(0), stationaryAcceptCount(0), stationaryCount(0), 
                    nucleotideTransitionRates(6.0, 0.0), oldNucleotideTransitionRates(6.0, 0.0){

    // fill a vector alpha with all ones and generate a random draw of stationaries that sum to 1
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<double> alpha(4, 2.0);
    Probability::Dirichlet::rv(&rng, alpha, currentStationary);
    oldStationary = currentStationary;

    // set the prior for stationary
    stationaryPrior = Probability::Dirichlet::lnPdf(alpha, currentStationary);

    // now set the nucleotideTransitionRates and the ratePrior
    ratePrior = 0;
    for(int nucIter = 0; nucIter < nucleotideTransitionRates.size(); nucIter++){
        nucleotideTransitionRates[nucIter] = Probability::Gamma::rv(&rng, 2, 5);
        ratePrior += Probability::Gamma::lnPdf(2, 5,nucleotideTransitionRates[nucIter]);
    }

    oldStationaryPrior = stationaryPrior;
    oldRatePrior = ratePrior;
    dirty();
}

// if the proposed rateMatrix update is proposed, figure out what kind of update was made (rate or stationary) and 
// copy the current parameters and priors.
void RateMatrix::accept() {

    // copy parameters
    oldNucleotideTransitionRates = nucleotideTransitionRates;
    oldStationary = currentStationary;

    // copy priors
    oldRatePrior = ratePrior;
    oldStationaryPrior = stationaryPrior;

    // now increment based on what type of update was made
    if(rateOrStationary == 2){
        rateAcceptCount++;
    }
    else if (rateOrStationary == 1){
        stationaryAcceptCount++;
    }
}

// if the proposed rateMatrix update is rejected, copy the old parameters and priors to current to discard the proposal
void RateMatrix::reject() {

    // copy parameters
    nucleotideTransitionRates = oldNucleotideTransitionRates;
    currentStationary = oldStationary;
    
    // copy priors
    ratePrior = oldRatePrior;
    stationaryPrior = oldStationaryPrior;

    // reset the update-type flag for tuning
    rateOrStationary = 0;
}

// returns the prior of the latest proposed parameters that have not been rejected yet
double RateMatrix::lnPrior() {
    return stationaryPrior + ratePrior;
}

// pick a single transition rate within the matrix and propose a rescaling update.
double RateMatrix::updateRates(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    rateCount++; // increment update type for tuning
    dirty(); // set flag for conditional likelihood recalculation
    rateOrStationary = 2; // set update-type flag to indicate a rate type update for tuning

    // randomly pick one transition rate to update
    int randIndex = Probability::Uniform::rv(&rng, 0, 6);
    double scale = std::exp(rateStepsize * (rng.uniformRv() - 0.5));
    nucleotideTransitionRates[randIndex] = nucleotideTransitionRates[randIndex] * scale;
    double log_hastings = std::log(scale);
    
    // now recalculate the scale prior
    ratePrior = 0;
    for(int nucIter = 0; nucIter < nucleotideTransitionRates.size(); nucIter++){
        ratePrior += Probability::Gamma::lnPdf(2, 5, nucleotideTransitionRates[nucIter]);
    }
    return log_hastings;
}

// make an update to the 4 nucleotide stationaries. 
double RateMatrix::updateStationary(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    dirty();
    stationaryCount++;
    rateOrStationary = 1; // set this flag to 1 to indicate that stationaries are being updated

    std::vector<double> z(currentStationary.size(), 0.0); // proposed new stationary vector
    std::vector<double> alphaForward(currentStationary.size(), 0.0); // vector for Dirichlet dist. params 
    std::vector<double> alphaReverse(currentStationary.size(), 0.0); // vector for Dirichlet dist. param but in reverse

    // calculate the alphaForward based on the current parameters.
    for(int i = 0; i < 4; i++){
        alphaForward[i] = (currentStationary[i] * stationaryAlpha) + 0.1;
    }

    // get the new stationary vector using alphaForward as the parameters
    Probability::Dirichlet::rv(&rng, alphaForward, z);

    // now calculate how to go backwards from 
    for(int i = 0; i < 4; i++){
        alphaReverse[i] = (z[i] * stationaryAlpha) + 0.1;
    }

    currentStationary = z; // save the new stationary vector
    double hastings = Probability::Dirichlet::lnPdf(alphaReverse, currentStationary) - Probability::Dirichlet::lnPdf(alphaForward, z);
    return hastings;
}

// This method calculates and returns the rateMatrix via using the stored nucleotideTransitionRates and stationary distribution
Matrix<double> RateMatrix::Q() {
    Matrix<double> returnMatrix(4, 4, 0.0);

    // fill the returnMatrix with the nucleotideTransitionRates first
    int index = 0;
    for(int i = 0; i < 4; i++){
        for(int j = i+1; j < 4; j++){
            returnMatrix(i,j) = nucleotideTransitionRates[index];
            returnMatrix(j,i) = nucleotideTransitionRates[index++];
        }
    }
    
    // now multiply the matrix by the stationary, make the diagonals the additive inverse
    double scaler= 0.0;
    for(int i = 0; i < 4; i++){
        double total = 0.0;
        for(int j = 0; j < 4; j++){

            // if not diagonal, multiply by stationary and add to total
            if(j!=i){
                returnMatrix(i,j) = currentStationary[j] * returnMatrix(i,j);
                total += returnMatrix(i , j);
            } 
        }

        // make diagonal elements the additive inverse for the row
        returnMatrix(i, i) = total * -1;
        scaler += returnMatrix(i, i) * currentStationary[i];
    }

    // scale matrix
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            returnMatrix(i, j) /= -1 * scaler;
    
    return returnMatrix;
}

// adjust the step-sizes for both rate and stationary updates based on acceptance rates during burn-in
void RateMatrix::tune(){ 
    double transitionRateRate = (double) rateAcceptCount/rateCount;
    if(transitionRateRate > 0.33){
        rateStepsize *= (1.0 + ((transitionRateRate-0.33)/0.67));
    }
    else {
        rateStepsize /= (2.0 - transitionRateRate/0.33);
    }

    // reset counts to evaluate next round of tuning
    rateAcceptCount = 0;
    rateCount = 0;
    std::cout << "transitionRate: " << transitionRateRate << "\n"; // print out for logging purposes 

    // we need to modify the stationaryAlpha
    double stationaryRate = (double) stationaryAcceptCount/stationaryCount;

    // if our acceptance rate is too high, that means we are exploring too much, so dial it back a little bit
    if(stationaryRate > 0.33){
        stationaryAlpha /= (1.0 + ((stationaryRate-0.33)/0.67));
    }
    else {
        stationaryAlpha *= (2.0 + ((stationaryRate-0.33)));
    }

    // reset counts to evaluate next round of tuning
    stationaryAcceptCount = 0;
    stationaryCount = 0;
    std::cout << "stationaryRate: " << stationaryRate << "\n"; // print out for logging purposes 
}