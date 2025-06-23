#ifndef RATE_EIGEN_HPP
#define RATE_EIGEN_HPP
#include <complex>
#include <vector>

struct RateEigen {
    double* c_ijk;
    double* eigenvalue;
    double* oldC_ijk;
    double* oldEigenvalue;
    
    // johndu addition to support gammaDist model
    double* diagLeftMatrix;
    double* oldDiagLeftMatrix;
    double* diagRightMatrix;
    double* oldDiagRightMatrix;
    int numStates;

    RateEigen(int nS) : numStates(nS) {
        c_ijk = new double[numStates * numStates * numStates];
        eigenvalue = new double[numStates];
        oldC_ijk = new double[numStates * numStates * numStates];
        oldEigenvalue = new double[numStates];
        diagLeftMatrix = new double[numStates * numStates];
        oldDiagLeftMatrix = new double[numStates * numStates];
        diagRightMatrix = new double[numStates * numStates];
        oldDiagRightMatrix = new double[numStates * numStates];

        std::fill(c_ijk, c_ijk + numStates * numStates * numStates, 0.0);
        std::fill(oldC_ijk, oldC_ijk + numStates * numStates * numStates, 0.0);
        std::fill(eigenvalue, eigenvalue + numStates, 0.0);
        std::fill(oldEigenvalue, oldEigenvalue + numStates, 0.0);
        std::fill(diagLeftMatrix, diagLeftMatrix + numStates * numStates, 0.0);
        std::fill(oldDiagLeftMatrix, oldDiagLeftMatrix + numStates * numStates, 0.0);
        std::fill(diagRightMatrix, diagRightMatrix + numStates * numStates, 0.0);
        std::fill(oldDiagRightMatrix, oldDiagRightMatrix + numStates * numStates, 0.0);
    }

    RateEigen(const RateEigen& other) {
        numStates = other.numStates;
        c_ijk = new double[numStates * numStates * numStates];
        eigenvalue = new double[numStates];
        oldC_ijk = new double[numStates * numStates * numStates];
        oldEigenvalue = new double[numStates];
        diagLeftMatrix = new double[numStates * numStates];
        oldDiagLeftMatrix = new double[numStates * numStates];
        diagRightMatrix = new double[numStates * numStates];
        oldDiagRightMatrix = new double[numStates * numStates];
        
        std::copy(other.c_ijk, other.c_ijk + numStates * numStates * numStates, c_ijk);
        std::copy(other.oldC_ijk, other.oldC_ijk + numStates * numStates * numStates, oldC_ijk);
        std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
        std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
        std::copy(other.diagLeftMatrix, other.diagLeftMatrix + numStates * numStates, diagLeftMatrix);
        std::copy(other.oldDiagLeftMatrix, other.oldDiagLeftMatrix + numStates * numStates, oldDiagLeftMatrix);
        std::copy(other.diagRightMatrix, other.diagRightMatrix + numStates * numStates, diagRightMatrix);
        std::copy(other.oldDiagRightMatrix, other.diagLeftMatrix + numStates * numStates, diagRightMatrix);
    }

    RateEigen& operator=(const RateEigen& other) {
        if (this != &other) {
            delete [] c_ijk;
            delete [] eigenvalue;
            delete [] oldC_ijk;
            delete [] oldEigenvalue;
            delete [] diagLeftMatrix;
            delete [] oldDiagLeftMatrix;
            delete [] oldDiagRightMatrix;
            delete [] oldDiagRightMatrix;


            numStates = other.numStates;
            c_ijk = new double[numStates * numStates * numStates];
            eigenvalue = new double[numStates];
            oldC_ijk = new double[numStates * numStates * numStates];
            oldEigenvalue = new double[numStates];
            diagLeftMatrix = new double[numStates * numStates];
            oldDiagLeftMatrix = new double[numStates * numStates];
            diagRightMatrix = new double[numStates * numStates];
            oldDiagRightMatrix = new double[numStates * numStates];

            std::copy(other.c_ijk, other.c_ijk + numStates * numStates * numStates, c_ijk);
            std::copy(other.oldC_ijk, other.oldC_ijk + numStates * numStates * numStates, oldC_ijk);
            std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
            std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
            std::copy(other.diagLeftMatrix, other.diagLeftMatrix + numStates * numStates, diagLeftMatrix);
            std::copy(other.oldDiagLeftMatrix, other.oldDiagLeftMatrix + numStates * numStates, oldDiagLeftMatrix);
            std::copy(other.diagRightMatrix, other.diagRightMatrix + numStates * numStates, diagRightMatrix);
            std::copy(other.oldDiagRightMatrix, other.diagLeftMatrix + numStates * numStates, diagRightMatrix);
        }
        return *this;
    }

    ~RateEigen() {
    delete [] c_ijk;
    delete [] eigenvalue;
    delete [] oldEigenvalue;
    delete [] oldC_ijk;
    delete [] diagLeftMatrix;
    delete [] oldDiagLeftMatrix;
    delete [] diagRightMatrix;
    delete [] oldDiagRightMatrix;
    }
};

struct ComplexRateEigen {
    std::complex<double>* cc_ijk; 
    std::complex<double>* ceigenvalue;
    std::complex<double>* oldCC_ijk;
    std::complex<double>* oldCeigenvalue;

    // johndu addition to support gammaDist model
    std::complex<double>* cDiagLeftMatrix;
    std::complex<double>* cDiagRightMatrix;
    std::complex<double>* oldCDiagLeftMatrix;
    std::complex<double>* oldCDiagRightMatrix;
    int numStates;

    ComplexRateEigen(int nS) : numStates(nS) {
        cc_ijk = new std::complex<double>[numStates * numStates * numStates];
        ceigenvalue = new std::complex<double>[numStates];
        oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        cDiagLeftMatrix = new std::complex<double>[numStates * numStates];
        cDiagRightMatrix = new std::complex<double>[numStates * numStates];
        oldCDiagLeftMatrix = new std::complex<double>[numStates * numStates];
        oldCDiagRightMatrix = new std::complex<double>[numStates * numStates];

        std::fill(cc_ijk, cc_ijk + numStates * numStates * numStates, 0.0);
        std::fill(oldCC_ijk, oldCC_ijk + numStates * numStates * numStates, 0.0);
        std::fill(ceigenvalue, ceigenvalue + numStates, 0.0);
        std::fill(oldCeigenvalue, oldCeigenvalue + numStates, 0.0);  
        std::fill(cDiagLeftMatrix, cDiagLeftMatrix + numStates * numStates, 0.0);
        std::fill(cDiagRightMatrix, cDiagRightMatrix + numStates * numStates, 0.0);
        std::fill(oldCDiagLeftMatrix, oldCDiagLeftMatrix + numStates * numStates, 0.0);
        std::fill(oldCDiagRightMatrix, oldCDiagRightMatrix + numStates * numStates, 0.0);
    }


    ComplexRateEigen(const ComplexRateEigen& other) {
        numStates = other.numStates;
        cc_ijk = new std::complex<double>[numStates * numStates * numStates];
        ceigenvalue = new std::complex<double>[numStates];
        oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        cDiagLeftMatrix = new std::complex<double>[numStates * numStates];
        cDiagRightMatrix = new std::complex<double>[numStates * numStates];
        oldCDiagLeftMatrix = new std::complex<double>[numStates * numStates];
        oldCDiagRightMatrix = new std::complex<double>[numStates * numStates];
        
        std::copy(other.cc_ijk, other.cc_ijk + numStates * numStates * numStates, cc_ijk);
        std::copy(other.oldCC_ijk, other.oldCC_ijk + numStates * numStates * numStates, oldCC_ijk);
        std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
        std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);
        std::copy(other.cDiagLeftMatrix, other.cDiagLeftMatrix + numStates * numStates, cDiagLeftMatrix);
        std::copy(other.cDiagRightMatrix, other.cDiagRightMatrix + numStates * numStates, cDiagRightMatrix);
        std::copy(other.oldCDiagLeftMatrix, other.oldCDiagLeftMatrix + numStates * numStates, oldCDiagLeftMatrix);
        std::copy(other.oldCDiagRightMatrix, other.oldCDiagRightMatrix + numStates * numStates, oldCDiagRightMatrix);
    }

    ComplexRateEigen& operator=(const ComplexRateEigen& other) {
        if (this != &other) {
            delete [] cc_ijk;
            delete [] ceigenvalue;
            delete [] oldCC_ijk;
            delete [] oldCeigenvalue;
            delete [] cDiagLeftMatrix;
            delete [] cDiagRightMatrix;
            delete [] oldCDiagLeftMatrix;
            delete [] oldCDiagRightMatrix;

            numStates = other.numStates;
            cc_ijk = new std::complex<double>[numStates * numStates * numStates];
            ceigenvalue = new std::complex<double>[numStates];
            oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
            oldCeigenvalue = new std::complex<double>[numStates];
            cDiagLeftMatrix = new std::complex<double>[numStates * numStates];
            cDiagRightMatrix = new std::complex<double>[numStates * numStates];
            oldCDiagLeftMatrix = new std::complex<double>[numStates * numStates];
            oldCDiagRightMatrix = new std::complex<double>[numStates * numStates];

            std::copy(other.cc_ijk, other.cc_ijk + numStates * numStates * numStates, cc_ijk);
            std::copy(other.oldCC_ijk, other.oldCC_ijk + numStates * numStates * numStates, oldCC_ijk);
            std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
            std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);
            std::copy(other.cDiagLeftMatrix, other.cDiagLeftMatrix + numStates * numStates, cDiagLeftMatrix);
            std::copy(other.cDiagRightMatrix, other.cDiagRightMatrix + numStates * numStates, cDiagRightMatrix);
            std::copy(other.oldCDiagLeftMatrix, other.oldCDiagLeftMatrix + numStates * numStates, oldCDiagLeftMatrix);
            std::copy(other.oldCDiagRightMatrix, other.oldCDiagRightMatrix + numStates * numStates, oldCDiagRightMatrix);
        }
        return *this;
    }

    ~ComplexRateEigen() {
    delete [] cc_ijk;
    delete [] ceigenvalue;
    delete [] oldCC_ijk;
    delete [] oldCeigenvalue;
    delete [] cDiagLeftMatrix;
    delete [] cDiagRightMatrix;
    delete [] oldCDiagLeftMatrix;
    delete [] oldCDiagRightMatrix;
    }
};

#endif