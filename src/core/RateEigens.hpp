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
    Matrix<double>* diagLeftMatrix;
    Matrix<double>* oldDiagLeftMatrix;
    Matrix<double>* diagRightMatrix;
    Matrix<double>* oldDiagRightMatrix;
    int numStates;

    RateEigen(int nS) : numStates(nS) {
        c_ijk = new double[numStates * numStates * numStates];
        eigenvalue = new double[numStates];
        oldC_ijk = new double[numStates * numStates * numStates];
        oldEigenvalue = new double[numStates];
        diagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        diagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);

        std::fill(c_ijk, c_ijk + numStates * numStates * numStates, 0.0);
        std::fill(oldC_ijk, oldC_ijk + numStates * numStates * numStates, 0.0);
        std::fill(eigenvalue, eigenvalue + numStates, 0.0);
        std::fill(oldEigenvalue, oldEigenvalue + numStates, 0.0);
    }

    RateEigen(const RateEigen& other) {
        numStates = other.numStates;
        c_ijk = new double[numStates * numStates * numStates];
        eigenvalue = new double[numStates];
        oldC_ijk = new double[numStates * numStates * numStates];
        oldEigenvalue = new double[numStates];
        diagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        diagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
        
        std::copy(other.c_ijk, other.c_ijk + numStates * numStates * numStates, c_ijk);
        std::copy(other.oldC_ijk, other.oldC_ijk + numStates * numStates * numStates, oldC_ijk);
        std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
        std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
        diagLeftMatrix->inject(*(other.diagLeftMatrix));
        oldDiagLeftMatrix->inject(*(other.oldDiagLeftMatrix));
        diagRightMatrix->inject(*(other.diagRightMatrix));
        oldDiagRightMatrix->inject(*(other.oldDiagRightMatrix));
    }

    RateEigen& operator=(const RateEigen& other) {
        if (this != &other) {
            delete [] c_ijk;
            delete [] eigenvalue;
            delete [] oldC_ijk;
            delete [] oldEigenvalue;
            delete diagLeftMatrix;
            delete diagRightMatrix;
            delete oldDiagLeftMatrix;
            delete oldDiagRightMatrix;

            numStates = other.numStates;
            c_ijk = new double[numStates * numStates * numStates];
            eigenvalue = new double[numStates];
            oldC_ijk = new double[numStates * numStates * numStates];
            oldEigenvalue = new double[numStates];

            diagLeftMatrix->inject(*(other.diagLeftMatrix));
            oldDiagLeftMatrix->inject(*(other.oldDiagLeftMatrix));
            diagRightMatrix->inject(*(other.diagRightMatrix));
            oldDiagRightMatrix->inject(*(other.oldDiagRightMatrix));

            std::copy(other.c_ijk, other.c_ijk + numStates * numStates * numStates, c_ijk);
            std::copy(other.oldC_ijk, other.oldC_ijk + numStates * numStates * numStates, oldC_ijk);
            std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
            std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
        }
        return *this;
    }

    ~RateEigen() {
    delete [] c_ijk;
    delete [] eigenvalue;
    delete [] oldEigenvalue;
    delete [] oldC_ijk;
    delete diagLeftMatrix;
    delete diagRightMatrix;
    delete oldDiagLeftMatrix;
    delete oldDiagRightMatrix;
    }
};

struct ComplexRateEigen {
    std::complex<double>* cc_ijk; 
    std::complex<double>* ceigenvalue;
    std::complex<double>* oldCC_ijk;
    std::complex<double>* oldCeigenvalue;

    // johndu addition to support gammaDist model
    Matrix<std::complex<double>>* cDiagLeftMatrix;
    Matrix<std::complex<double>>* cDiagRightMatrix;
    Matrix<std::complex<double>>* oldCDiagLeftMatrix;
    Matrix<std::complex<double>>* oldCDiagRightMatrix;
    int numStates;

    ComplexRateEigen(int nS) : numStates(nS) {
        cc_ijk = new std::complex<double>[numStates * numStates * numStates];
        ceigenvalue = new std::complex<double>[numStates];
        oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        cDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        cDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);

        std::fill(cc_ijk, cc_ijk + numStates * numStates * numStates, 0.0);
        std::fill(oldCC_ijk, oldCC_ijk + numStates * numStates * numStates, 0.0);
        std::fill(ceigenvalue, ceigenvalue + numStates, 0.0);
        std::fill(oldCeigenvalue, oldCeigenvalue + numStates, 0.0); 
    }


    ComplexRateEigen(const ComplexRateEigen& other) {
        numStates = other.numStates;
        cc_ijk = new std::complex<double>[numStates * numStates * numStates];
        ceigenvalue = new std::complex<double>[numStates];
        oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        cDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        cDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        
        std::copy(other.cc_ijk, other.cc_ijk + numStates * numStates * numStates, cc_ijk);
        std::copy(other.oldCC_ijk, other.oldCC_ijk + numStates * numStates * numStates, oldCC_ijk);
        std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
        std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);

        cDiagLeftMatrix->inject(*(other.cDiagLeftMatrix));
        oldCDiagLeftMatrix->inject(*(other.oldCDiagLeftMatrix));
        cDiagRightMatrix->inject(*(other.cDiagRightMatrix));
        oldCDiagRightMatrix->inject(*(other.oldCDiagRightMatrix));
    }

    ComplexRateEigen& operator=(const ComplexRateEigen& other) {
        if (this != &other) {
            delete [] cc_ijk;
            delete [] ceigenvalue;
            delete [] oldCC_ijk;
            delete [] oldCeigenvalue;
            delete cDiagLeftMatrix;
            delete cDiagRightMatrix;
            delete oldCDiagLeftMatrix;
            delete oldCDiagRightMatrix;

            numStates = other.numStates;
            cc_ijk = new std::complex<double>[numStates * numStates * numStates];
            ceigenvalue = new std::complex<double>[numStates];
            oldCC_ijk = new std::complex<double>[numStates * numStates * numStates];
            oldCeigenvalue = new std::complex<double>[numStates];
            cDiagLeftMatrix->inject(*(other.cDiagLeftMatrix));
            cDiagRightMatrix->inject(*(other.cDiagRightMatrix));
            oldCDiagLeftMatrix->inject(*(other.oldCDiagLeftMatrix));
            oldCDiagRightMatrix->inject(*(other.oldCDiagRightMatrix));

            std::copy(other.cc_ijk, other.cc_ijk + numStates * numStates * numStates, cc_ijk);
            std::copy(other.oldCC_ijk, other.oldCC_ijk + numStates * numStates * numStates, oldCC_ijk);
            std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
            std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);
        }
        return *this;
    }

    ~ComplexRateEigen() {
    delete [] cc_ijk;
    delete [] ceigenvalue;
    delete [] oldCC_ijk;
    delete [] oldCeigenvalue;
    delete cDiagLeftMatrix;
    delete cDiagRightMatrix;
    delete oldCDiagLeftMatrix;
    delete oldCDiagRightMatrix;
    }
};

#endif