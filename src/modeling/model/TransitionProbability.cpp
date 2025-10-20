#include <complex>
#include "TransitionProbability.hpp"
#include "core/RateEigens.hpp"
#include "core/Math.hpp"
#include <cstring>
#include "test.h"

TransitionProbability::TransitionProbability(const int nn)
    : numStates(4), numNodes(nn), probs1(), probs2() {

	Matrix<double> Q(numStates, numStates, 0.0);
	eigens = new EigenSystem(numStates);
	
	allocateQ(1);
	updateQ(Q);
}

/* Destructor. Deallocates memory used for Q matrix
   and eigensystem. */
TransitionProbability::~TransitionProbability(void) {
	delete eigens;
	
	for(auto i : probs1){
		delete [] i; 
	}
	for(auto i : probs2){
		delete [] i;
	}
}

/* This method calculates the transition probability for the branch that flows into param node. Stores the calculated transition
   probability matrix in a buffer based on param state, rate and node. Rate is an artifact of a legacy model. State indicates which
   buffer to use in the ping-pong buffer and each node has its own ancestor branch and therefore has a transition probabiliy associated
   with the node. Param alpha and beta are parameters of the gamma distribution that governs the branch length of the ancestor branch.
*/
void TransitionProbability::setProbs(const int state, const int rate, const int node, double alpha, double beta) {

	// fetch the matrix based off params, state, rate and node. This matrix will store transition probability
	Matrix<double> P0 = (*this)(state, rate, node);  
 
	// call the method to transform P0
	tiProbsGamma(alpha, beta, P0);
}

// Uses formula outlined in Huelsenbeck's "Bayesian Perspective on a Non-parsimonious Parsimony Model" to
// calculate transition probability matrix ((IdentityMatrix - (1/rate|beta) * rateMatrix)^-(shape|alpha))
// using the rate matrix, Q, which is stored as a class static variable.
void TransitionProbability::tiProbsGamma(const double shape, const double rate, Matrix<double>& P0) {

	// first copy the rate matrix so we can transform it
	Matrix<double> transformMatrix(Q.copy());

	// next perform a scalar multiplication by 1/rate and then subtract the transformed matrix from the identity matrix (I-A)
	transformMatrix *= -1/rate;
	for(int i = 0; i < transformMatrix.dim1(); i++){
		transformMatrix(i,i) = 1 + transformMatrix(i,i);
	}

	// get the eigenvalues and eigenvectors via eigens->update. x
	isComplex = eigens->update(transformMatrix, rateEigen[0], complexRateEigen[0]);
	
	// now split here depending on whether or not we have complex eigenvalues
	if(!isComplex){
		RateEigen newDiag = rateEigen[0];
		Matrix<double> *leftMatrix = newDiag.diagLeftMatrix;
		Matrix<double> *rightMatrix = newDiag.diagRightMatrix;
		double* eigenvalues = newDiag.eigenvalue;

		// take the eigenvalues to the power of -(shape|alpha) to mimic A^-x = P * D^-x * P^-1 as part of taking the matrix to the
		// -(shape|alpha) power to get our diagonal matrix
		Matrix<double> diagonalMatrix(Q.dim1(), Q.dim2(), 0.0);
		for(int i = 0; i < Q.dim1(); i++){
			diagonalMatrix(i, i) = eigenvalues[i];
		}

		// take the diagonal matrix to the -shape power. Here we use the property a^b = e^(b * log(a)) 
		for(int i = 0; i < Q.dim1(); i++){
			diagonalMatrix(i, i) = std::pow(eigenvalues[i], -1 * shape);
		}

		// now get transitionProbability matrix by using the property A = P*D*P^-1
		Matrix<double> newMatrix = ((*leftMatrix) * diagonalMatrix);
		Matrix<double> finalMatrix = newMatrix * (*rightMatrix);
		
		// copy over the matrix to designated buffer for transition probability matrix
		for(int i = 0; i < Q.dim1(); i++){
			for(int j = 0; j < Q.dim2(); j++){
				P0(i,j) = finalMatrix(i, j); 
			}
		}
	} else {
		ComplexRateEigen newDiag = complexRateEigen[0];
		Matrix<std::complex<double>> *rightMatrix = newDiag.cDiagRightMatrix;
		Matrix<std::complex<double>> *leftMatrix = newDiag.cDiagLeftMatrix;
		std::complex<double> *cEigenvalues = newDiag.ceigenvalue;
		
		// now take the complex numbers to the power of -1/shape
		Matrix<std::complex<double>> diagonalMatrix(Q.dim1(), Q.dim2(), 0.0);
		for(int i = 0; i < Q.dim1(); i++){
			diagonalMatrix(i, i) = std::pow(cEigenvalues[i], -1 * shape);
		}
		
		// now that the transition probability matrix has been calculated, because a matrix of real numbers ^ real positive numbers
		// is still a matrix of real numbers, so we can safely discard the imaginary components of the transition probability
		Matrix<std::complex<double>> newMatrix = (( (*leftMatrix) * diagonalMatrix)*(*rightMatrix));  
		for(int i = 0; i < Q.dim1(); i++){
			for(int j = 0; j < Q.dim2(); j++){
				P0(i,j) = newMatrix(i, j).real(); 
			}
		}
	} 
}

// a simple method to update the rateMatrix Q
void TransitionProbability::updateQ(Matrix<double> otherQ){
	Q = otherQ;
}

void TransitionProbability::allocateQ(int size){

	rateEigen.push_back(RateEigen(numStates));
	complexRateEigen.push_back(ComplexRateEigen(numStates)); //initialize a ComplexRateEigen struct

	probs1.push_back(new Matrix<double>[numNodes]);
	probs2.push_back(new Matrix<double>[numNodes]);

	// for each node, initialize two matrices, one in probs1 and probs 2
	for(int j = 0; j < numNodes; j++){
		probs1.back()[j] = Matrix<double>(numStates, numStates, 0.0);
		probs2.back()[j] = Matrix<double>(numStates, numStates, 0.0);
	}

	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}

// Be sure you want to delete!!
void TransitionProbability::deleteQ(const int index) {
	complexRateEigen.erase(complexRateEigen.begin() + index);

	auto prob_it1 = probs1.begin() + index;
	delete [] *prob_it1;
	probs1.erase(prob_it1);

	auto prob_it2 = probs2.begin() + index;
	delete [] *prob_it2;
	probs2.erase(prob_it2);

	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}

void TransitionProbability::deleteNQ(const int count) {
	for(int i = 0; i < count; i++){
		rateEigen.pop_back();
		complexRateEigen.pop_back();

		auto probs_it1 = std::prev(probs1.end());
		delete [] *probs_it1;
		probs1.pop_back();

		auto probs_it2 = std::prev(probs2.end());
		delete [] *probs_it2;
		probs2.pop_back();
	}
	
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}