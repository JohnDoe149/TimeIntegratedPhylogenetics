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
	accept();
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

/* Store the diagonalized matrix via its eigenvalues and eigenvectors. Move the current eigenvalues and eigenvectors to
   the "old" version to store them in case a backwards move is necessary.
*/
void TransitionProbability::accept(void) {
	isOldComplex = isComplex;
	for(int i = 0; i < isComplex.size(); i++){
		if(!isComplex[i]){
			memcpy(rateEigen[i].oldEigenvalue, rateEigen[i].eigenvalue, numStates*sizeof(double));
            rateEigen[i].oldDiagLeftMatrix->inject(*rateEigen[i].diagLeftMatrix);
            rateEigen[i].oldDiagRightMatrix->inject(*rateEigen[i].diagRightMatrix);
		}
		else {
			memcpy(complexRateEigen[i].oldCeigenvalue, complexRateEigen[i].ceigenvalue, numStates*sizeof(std::complex<double>));
            complexRateEigen[i].oldCDiagLeftMatrix->inject(*complexRateEigen[i].cDiagLeftMatrix);
            complexRateEigen[i].oldCDiagRightMatrix->inject(*complexRateEigen[i].cDiagRightMatrix);
		}
	}
}

/* If the proposed transition probability is rejected, then the previous accepted transiiton probability must be reused.  
   In other words, copy the eigenvalues and eigenvectors that make up the diagonalized matrix from the "old" or "storage" versions
   of the variables into the current versions of the variables.
*/
void TransitionProbability::reject(void) {	
	isComplex = isOldComplex;
	for(int i = 0; i < isOldComplex.size(); i++){
		if(!isComplex[i]){
			memcpy(rateEigen[i].eigenvalue, rateEigen[i].oldEigenvalue, numStates*sizeof(double));
			memcpy(rateEigen[i].diagLeftMatrix, rateEigen[i].oldDiagLeftMatrix, numStates*numStates*sizeof(double));
			memcpy(rateEigen[i].diagRightMatrix, rateEigen[i].oldDiagRightMatrix, numStates*numStates*sizeof(double));
		}
		else {
			memcpy(complexRateEigen[i].ceigenvalue, complexRateEigen[i].oldCeigenvalue, numStates*sizeof(std::complex<double>));
			memcpy(complexRateEigen[i].cDiagLeftMatrix, complexRateEigen[i].oldCDiagLeftMatrix, numStates*numStates*sizeof(std::complex<double>));
			memcpy(complexRateEigen[i].cDiagRightMatrix, complexRateEigen[i].oldCDiagRightMatrix, numStates*numStates*sizeof(std::complex<double>));
		}
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
// calculate transition probability matrix ((IdentityMatrix - (1/scale|beta) * rateMatrix)^-(shape|alpha))
// using the rate matrix, Q, which is stored as a class static variable.
void TransitionProbability::tiProbsGamma(const double shape, const double scale, Matrix<double>& P0) {

	// first copy the rate matrix so we can transform it
	Matrix<double> transformMatrix(Q.copy());

	#ifdef TRANSPROB_PRINT
	std::cout << "starting matrix: \n";
	transformMatrix.print();
	#endif

	// next perform a scalar multiplication by 1/scale and then subtract the transformed matrix from the identity matrix (I-A)
	transformMatrix *= -1/scale;
	for(int i = 0; i < transformMatrix.dim1(); i++){
		transformMatrix(i,i) = 1 + transformMatrix(i,i);
	}

	#ifdef TRANSPROB_PRINT
	std::cout << "matrix predecomposition: \n";
	transformMatrix.print();
	#endif

	// get the eigenvalues and eigenvectors via eigens->update. x
	isComplex[0] = eigens->update(transformMatrix, rateEigen[0], complexRateEigen[0]);
	
	// now split here depending on whether or not we have complex eigenvalues
	if(!isComplex[0]){
		RateEigen newDiag = rateEigen[0];
		Matrix<double> *leftMatrix = newDiag.diagLeftMatrix;
		Matrix<double> *rightMatrix = newDiag.diagRightMatrix;
		double* eigenvalues = newDiag.eigenvalue;

		#ifdef TRANSPROB_PRINT
		std::cout << "rightmatrix postdecomposition: \n";
		rightMatrix->print();
		std::cout << "left matrix postdecomposition: \n";
		leftMatrix->print();
		std::cout << "eigenvalues: \n";
		std::cout << eigenvalues[0] << " " << eigenvalues[1] << " " << eigenvalues[2] << " " << eigenvalues[3] << "\n" << std::flush;
		#endif

		// take the eigenvalues to the power of -(shape|alpha) to mimic A^-x = P * D^-x * P^-1 as part of taking the matrix to the
		// -(shape|alpha) power to get our diagonal matrix
		Matrix<double> diagonalMatrix(Q.dim1(), Q.dim2(), 0.0);
		for(int i = 0; i < Q.dim1(); i++){
			diagonalMatrix(i, i) = eigenvalues[i];
		}

		#ifdef TRANSPROB_PRINT
		std::cout << "diagonal matrix \n";
		diagonalMatrix.print();
		std::cout << "verify correct diagonalization \n";
		Matrix<double> verifyMatrix = (((*leftMatrix) * diagonalMatrix)*(*rightMatrix));  
		verifyMatrix.print();
		#endif

		// take the diagonal matrix to the -shape power. Here we use the property a^b = e^(b * log(a)) 
		for(int i = 0; i < Q.dim1(); i++){
			diagonalMatrix(i, i) = std::exp(-shape * std::log(eigenvalues[i]));
		}

		#ifdef TRANSPROB_PRINT
		std::cout << "diagonal matrix taken to the power: \n";
		diagonalMatrix.print();
		#endif

		// now get transitionProbability matrix by using the property A = P*D*P^-1
		Matrix<double> newMatrix = ((*leftMatrix) * diagonalMatrix);
		Matrix<double> finalMatrix = newMatrix * (*rightMatrix);

		#ifdef TRANSPROB_PRINT
		std::cout << "P * D\n";
		newMatrix.print();
		std::cout << "PD * P^-1\n";
		finalMatrix.print();
		#endif
		
		// copy over the matrix to designated buffer for transition probability matrix
		for(int i = 0; i < Q.dim1(); i++){
			for(int j = 0; j < Q.dim2(); j++){
				P0(i,j) = finalMatrix(i, j); 
			}
		}

		#ifdef TRANSPROB_PRINT
		std::cout << "transition probability matrix: \n";
		P0.print();
		#endif
	} else {
		ComplexRateEigen newDiag = complexRateEigen[0];
		Matrix<std::complex<double>> *rightMatrix = newDiag.cDiagRightMatrix;
		Matrix<std::complex<double>> *leftMatrix = newDiag.cDiagLeftMatrix;
		std::complex<double> *cEigenvalues = newDiag.ceigenvalue;

		#ifdef TRANSPROB_PRINT
		std::cout << "complex rightmatrix postdecomposition: \n";
		rightMatrix->print();
		std::cout << "complex left matrix postdecomposition: \n";
		leftMatrix->print();		
		std::cout << "complex eigenvalues: \n";
		std::cout << cEigenvalues[0].real() << "+" << cEigenvalues[0].imag() << "i " << cEigenvalues[1].real() << "+" << cEigenvalues[1].imag() << "i " 
				  << cEigenvalues[2].real() << "+" << cEigenvalues[2].imag() << "i " << cEigenvalues[3].real() << "+" << cEigenvalues[3].imag() << "i "<< "\n" << std::flush;
		#endif
		
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

		#ifndef TRANSPROB_PRINT
		std::cout << "transition probability matrix: \n";
		P0.print();
		#endif
	} 
}

// a simple method to update the rateMatrix Q
void TransitionProbability::updateQ(Matrix<double> otherQ){
	Q = otherQ;
}

void TransitionProbability::allocateQ(int size){
	if(size > isComplex.size()) {
		for(int i = 0, num = size - isComplex.size(); i < num; i++){
			isComplex.push_back(false);
			rateEigen.push_back(RateEigen(numStates));
			complexRateEigen.push_back(ComplexRateEigen(numStates)); //initialize a ComplexRateEigen struct

			probs1.push_back(new Matrix<double>[numNodes]);
			probs2.push_back(new Matrix<double>[numNodes]);

			// for each node, initialize two matrices, one in probs1 and probs 2
			for(int j = 0; j < numNodes; j++){
				probs1.back()[j] = Matrix<double>(numStates, numStates, 0.0);
       	 		probs2.back()[j] = Matrix<double>(numStates, numStates, 0.0);
			}
		}
	}

	isComplex.shrink_to_fit();
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}

// Be sure you want to delete!!
void TransitionProbability::deleteQ(const int index) {
	isComplex.erase(isComplex.begin() + index);
	rateEigen.erase(rateEigen.begin() + index);
	complexRateEigen.erase(complexRateEigen.begin() + index);

	auto prob_it1 = probs1.begin() + index;
	delete [] *prob_it1;
	probs1.erase(prob_it1);

	auto prob_it2 = probs2.begin() + index;
	delete [] *prob_it2;
	probs2.erase(prob_it2);

	isComplex.shrink_to_fit();
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}

void TransitionProbability::deleteNQ(const int count) {
	for(int i = 0; i < count; i++){
		isComplex.pop_back();
		rateEigen.pop_back();
		complexRateEigen.pop_back();

		auto probs_it1 = std::prev(probs1.end());
		delete [] *probs_it1;
		probs1.pop_back();

		auto probs_it2 = std::prev(probs2.end());
		delete [] *probs_it2;
		probs2.pop_back();
	}
	
	isComplex.shrink_to_fit();
	rateEigen.shrink_to_fit();
	complexRateEigen.shrink_to_fit();
	probs1.shrink_to_fit();
	probs2.shrink_to_fit();
}