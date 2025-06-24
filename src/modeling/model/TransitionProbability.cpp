#include <complex>
#include "TransitionProbability.hpp"
#include "core/RateEigens.hpp"
#include "core/Math.hpp"
#include <cstring>

TransitionProbability::TransitionProbability(const int nn)
    : numStates(4), numNodes(nn), probs1(), probs2() {

	/*
	probs[0] = new Matrix<double>*[2*numNodes*numCats];
    probs[1] = probs[0] + numNodes;

    for(int i = 0; i < numNodes*numCats; i++){
        probs[0][i] = new Matrix<double>(numStates, numStates, 0.0);
        probs[1][i] = new Matrix<double>(numStates, numStates, 0.0);
    }
	*/

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

void TransitionProbability::accept(void) {
	isOldComplex = isComplex;
	for(int i = 0; i < isComplex.size(); i++){
		if(!isComplex[i]){
			memcpy(rateEigen[i].oldC_ijk, rateEigen[i].c_ijk, numStates*numStates*numStates*sizeof(double));
			memcpy(rateEigen[i].oldEigenvalue, rateEigen[i].eigenvalue, numStates*sizeof(double));
		}
		else {
			memcpy(complexRateEigen[i].oldCC_ijk, complexRateEigen[i].cc_ijk, numStates*numStates*numStates*sizeof(std::complex<double>));
			memcpy(complexRateEigen[i].oldCeigenvalue, complexRateEigen[i].ceigenvalue, numStates*sizeof(std::complex<double>));
		}
	}
}

void TransitionProbability::reject(void) {	
	isComplex = isOldComplex;
	for(int i = 0; i < isOldComplex.size(); i++){
		if(!isComplex[i]){
			memcpy(rateEigen[i].c_ijk, rateEigen[i].oldC_ijk, numStates*numStates*numStates*sizeof(double));
			memcpy(rateEigen[i].eigenvalue, rateEigen[i].oldEigenvalue, numStates*sizeof(double));
		}
		else {
			memcpy(complexRateEigen[i].cc_ijk, complexRateEigen[i].oldCC_ijk, numStates*numStates*numStates*sizeof(std::complex<double>));
			memcpy(complexRateEigen[i].ceigenvalue, complexRateEigen[i].oldCeigenvalue, numStates*sizeof(std::complex<double>));
		}
	}
}

// FIX THIS TO CALCULATE TRANSITION PROB NEW WAY
void TransitionProbability::setProbs(const int state, const int rate, const int node, double alpha, double beta) {

	// state decides which buffer the matrix will be pulled from. Node and rate help index into the buffer
	// to get the right matrix. Rate is always 0, however due to the nature of how this model functions
	Matrix<double> P0 = (*this)(state, rate, node);  // this is the transition probability
	tiProbsGamma(alpha, beta, P0);
}

// Uses formula outlined in Huelsenbeck's "Bayesian Perspective on a Non-parsimonious Parsimony Model" to
// calculate transition probability matrix ((IdentityMatrix - (1/scale|beta) * rateMatrix)^-(shape|alpha))
void TransitionProbability::tiProbsGamma(const double shape, const double scale, Matrix<double>& P0) {
	Matrix<double> temp(Q.copy());
	std::cout << "starting matrix: \n";
	temp.print();
	temp *= 1/scale;
	for(int i = 0; i < temp.dim1(); i++){
		temp(i,i) = 1 - temp(i,i);
	}
	std::cout << "matrix predecomposition: \n";
	temp.print();
	// get the eigenvalues and eigenvectors
	Matrix<double> eigenDecompTemp(temp.copy());
	isComplex[0] = eigens->update(eigenDecompTemp, rateEigen[0], complexRateEigen[0]);
	RateEigen newEigenVectors = rateEigen[0];
	Matrix<double> *leftMatrix = newEigenVectors.diagLeftMatrix;
	Matrix<double> *rightMatrix = newEigenVectors.diagRightMatrix;
	std::cout << "rightmatrix postdecomposition: \n";
	rightMatrix->print();
	std::cout << "left matrix postdecomposition: \n";
	leftMatrix->print();
	double* eigenvalues = newEigenVectors.eigenvalue;
	std::cout << "eigenvalues: \n";
	std::cout << eigenvalues[0] << " " << eigenvalues[1] << " " << eigenvalues[2] << " " << eigenvalues[3] << "\n" << std::flush;
	// now print rateEigen and complexRateEigen and verify it worked correctly
	// also modify eigens->update to insert the eigenvectors into leftDiagMatrix
}

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