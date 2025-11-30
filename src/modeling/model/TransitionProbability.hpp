#ifndef TRANSITION_PROBABILITY_HPP
#define TRANSITION_PROBABILITY_HPP

#include <complex>
#include <vector>
#include "core/EigenSystem.hpp"
#include "core/Matrix.hpp"
#include "core/RateEigens.hpp"

class TransitionProbability {

	public:
                                TransitionProbability(const int nn);
                               ~TransitionProbability ();
        const Matrix<double>&   operator()(int s, int r, int n) const {
                                    return (s == 0) ? probs1[r][n] : probs2[r][n];
                                }
        int                     getNumStates(void) { return numStates; }          
        void                    setProbs(const int state, const int r, const int node, double alpha, double beta);
        std::vector<Matrix<double>>     generateProbs(Matrix<double> Q, std::vector<double> branches);
        void                    updateQ(Matrix<double> Q);
        void                    deleteQ(const int index);
        void                    deleteNQ(const int count);
        void                    allocateQ(int size);

    private:
        EigenSystem*            eigens; // eigens used to diagonalize matrices for matrix exponentiation
        std::vector<RateEigen>  rateEigen; // rateEigen stores a matrix's eigendecomposition for real values
        std::vector<ComplexRateEigen> complexRateEigen; // complexRateEigen stores a matrix's eigendecomposition for real values
        bool                    isComplex; // indicates whether the current eigensystem has complex eigenvalues
        int                     numNodes; // number of nodes in the tree
        int                     numStates; // number of states in the model (in this case 4, for nucleotides)
        Matrix<double>          Q; // Q matrix pulled from rate matrix
        std::vector<Matrix<double>*> probs1; 
        std::vector<Matrix<double>*> probs2;
        void                    tiProbsGamma(const double shape, const double scale, Matrix<double> &rateMatrix);
        void                    tiProbsComplexEigens(const double v, Matrix<double> &P, ComplexRateEigen& rE);
        void                    tiProbsEigens(const double v, Matrix<double> &P, RateEigen& rE);
};

#endif
