#ifndef RATE_EIGEN_HPP
#define RATE_EIGEN_HPP
#include <complex>
#include <vector>

// stores a matrix's eigendecomposition for real values
struct RateEigen {
    double* eigenvalue;
    double* oldEigenvalue;
    Matrix<double>* diagLeftMatrix;
    Matrix<double>* oldDiagLeftMatrix;
    Matrix<double>* diagRightMatrix;
    Matrix<double>* oldDiagRightMatrix;
    int numStates;

    RateEigen(int nS) : numStates(nS) {
        eigenvalue = new double[numStates];
        oldEigenvalue = new double[numStates];
        diagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        diagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);

        std::fill(eigenvalue, eigenvalue + numStates, 0.0);
        std::fill(oldEigenvalue, oldEigenvalue + numStates, 0.0);
    }

    RateEigen(const RateEigen& other) {
        numStates = other.numStates;
        eigenvalue = new double[numStates];
        oldEigenvalue = new double[numStates];
        diagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
        diagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
        oldDiagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
        
        std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
        std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
        diagLeftMatrix->inject(*(other.diagLeftMatrix));
        oldDiagLeftMatrix->inject(*(other.oldDiagLeftMatrix));
        diagRightMatrix->inject(*(other.diagRightMatrix));
        oldDiagRightMatrix->inject(*(other.oldDiagRightMatrix));
    }

    RateEigen& operator=(const RateEigen& other) {
        if (this != &other) {
            delete [] eigenvalue;
            delete [] oldEigenvalue;
            delete diagLeftMatrix;
            delete diagRightMatrix;
            delete oldDiagLeftMatrix;
            delete oldDiagRightMatrix;

            numStates = other.numStates;
            eigenvalue = new double[numStates];
            oldEigenvalue = new double[numStates];
            diagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
            oldDiagLeftMatrix = new Matrix<double>(numStates, numStates, 0.0);
            diagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);
            oldDiagRightMatrix = new Matrix<double>(numStates, numStates, 0.0);

            diagLeftMatrix->inject(*(other.diagLeftMatrix));
            oldDiagLeftMatrix->inject(*(other.oldDiagLeftMatrix));
            diagRightMatrix->inject(*(other.diagRightMatrix));
            oldDiagRightMatrix->inject(*(other.oldDiagRightMatrix));

            std::copy(other.eigenvalue, other.eigenvalue + numStates, eigenvalue);
            std::copy(other.oldEigenvalue, other.oldEigenvalue + numStates, oldEigenvalue);
        }
        return *this;
    }

    ~RateEigen() {
    delete [] eigenvalue;
    delete [] oldEigenvalue;
    delete diagLeftMatrix;
    delete diagRightMatrix;
    delete oldDiagLeftMatrix;
    delete oldDiagRightMatrix;
    }
};

// stores a matrix's eigendecomposition for complex values
struct ComplexRateEigen {
    std::complex<double>* ceigenvalue;
    std::complex<double>* oldCeigenvalue;

    // johndu addition to support gammaDist model
    Matrix<std::complex<double>>* cDiagLeftMatrix;
    Matrix<std::complex<double>>* cDiagRightMatrix;
    Matrix<std::complex<double>>* oldCDiagLeftMatrix;
    Matrix<std::complex<double>>* oldCDiagRightMatrix;
    int numStates;

    ComplexRateEigen(int nS) : numStates(nS) {
        ceigenvalue = new std::complex<double>[numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        cDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        cDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);

        std::fill(ceigenvalue, ceigenvalue + numStates, 0.0);
        std::fill(oldCeigenvalue, oldCeigenvalue + numStates, 0.0); 
    }


    ComplexRateEigen(const ComplexRateEigen& other) {
        numStates = other.numStates;
        ceigenvalue = new std::complex<double>[numStates];
        oldCeigenvalue = new std::complex<double>[numStates];
        cDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        cDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagLeftMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        oldCDiagRightMatrix = new Matrix<std::complex<double>>(numStates, numStates, 0.0);
        
        std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
        std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);

        cDiagLeftMatrix->inject(*(other.cDiagLeftMatrix));
        oldCDiagLeftMatrix->inject(*(other.oldCDiagLeftMatrix));
        cDiagRightMatrix->inject(*(other.cDiagRightMatrix));
        oldCDiagRightMatrix->inject(*(other.oldCDiagRightMatrix));
    }

    ComplexRateEigen& operator=(const ComplexRateEigen& other) {
        if (this != &other) {
            delete [] ceigenvalue;
            delete [] oldCeigenvalue;
            delete cDiagLeftMatrix;
            delete cDiagRightMatrix;
            delete oldCDiagLeftMatrix;
            delete oldCDiagRightMatrix;

            numStates = other.numStates;
            ceigenvalue = new std::complex<double>[numStates];
            oldCeigenvalue = new std::complex<double>[numStates];


            cDiagLeftMatrix->inject(*(other.cDiagLeftMatrix));
            cDiagRightMatrix->inject(*(other.cDiagRightMatrix));
            oldCDiagLeftMatrix->inject(*(other.oldCDiagLeftMatrix));
            oldCDiagRightMatrix->inject(*(other.oldCDiagRightMatrix));

            std::copy(other.ceigenvalue, other.ceigenvalue + numStates, ceigenvalue);
            std::copy(other.oldCeigenvalue, other.oldCeigenvalue + numStates, oldCeigenvalue);
        }
        return *this;
    }

    ~ComplexRateEigen() {
    delete [] ceigenvalue;
    delete [] oldCeigenvalue;
    delete cDiagLeftMatrix;
    delete cDiagRightMatrix;
    delete oldCDiagLeftMatrix;
    delete oldCDiagRightMatrix;
    }
};

#endif