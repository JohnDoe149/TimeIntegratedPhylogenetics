#include <gtest/gtest.h>
#include "core/Matrix.hpp"
#include "core/EigenSystem.hpp"
#include "core/RateEigens.hpp"


// test scalar multiplication
TEST(MatrixTest, Multiplication) {
  Matrix<double> testMatrix(4,4,1);
  Matrix<double> assertMatrix(4,4,1);

  testMatrix *= 0.137;
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      assertMatrix(i,j) = assertMatrix(i,j) * 0.137;
    }
  }

  EXPECT_TRUE(testMatrix == assertMatrix);
}

// verify that eigensystem diagonalization works correctly
TEST(MatrixTest, DiagonalTest){
  Matrix<double> testMatrix(4,4,1);
  
  // an arbitrary rate matrix to diagonalize nicely
  std::vector<double> stationary{0.25, 0.10, 0.35, 0.30};
  for(int i = 0; i < 4; i++){
    testMatrix(0, i) *= stationary[i];
    testMatrix(1, i) *= stationary[i];
    testMatrix(2, i) *= stationary[i];
    testMatrix(3, i) *= stationary[i];
  }
  for(int i = 0; i < 4; i++){
    double additiveInverse = 0;
    for(int j = 0; j < 4; j++){
      if(i!=j){
        additiveInverse += testMatrix(i,j);
      }
    }
    testMatrix(i,i) = -1 * additiveInverse;
  }

  // call diagonalize
  EigenSystem eigenSystem(4);
  RateEigen rateEigen(4);
  ComplexRateEigen complexRateEigen(4);
  eigenSystem.update(testMatrix, rateEigen, complexRateEigen);

  // now rateEigen should contain our diagonalized form, so check to
  Matrix<double> identityMatrix(4, 4, 0);
  for(int i = 0; i < 4; i++){
    identityMatrix(i,i) = rateEigen.eigenvalue[i];
  }

	Matrix<double> newMatrix = ((*rateEigen.diagLeftMatrix) * identityMatrix);
	Matrix<double> finalMatrix = newMatrix * (*rateEigen.diagRightMatrix);

  // some floating precision error so use EXPECT_NEAR
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      EXPECT_NEAR(finalMatrix(i,j), testMatrix(i,j), 0.0001);
    }
  }
}

// verify complex eigendecomposition works as expected
TEST(MatrixTest, ComplexDiagonalTest){
  Matrix<double> testMatrix(4,4,0.0);
  testMatrix(0,1) = -1;
  testMatrix(1,0) = 1;
  testMatrix(2,2) = 2;
  testMatrix(2,3) = -3;
  testMatrix(3,2) = 3;
  testMatrix(3,3) = 2;

  // call diagonalize
  EigenSystem eigenSystem(4);
  RateEigen rateEigen(4);
  ComplexRateEigen complexRateEigen(4);
  eigenSystem.update(testMatrix, rateEigen, complexRateEigen);

  Matrix<std::complex<double>> identityMatrix(4, 4, 0);
  for(int i = 0; i < 4; i++){
    identityMatrix(i,i) = complexRateEigen.ceigenvalue[i];
  }

	Matrix<std::complex<double>> newMatrix = ((*complexRateEigen.cDiagLeftMatrix) * identityMatrix);
	Matrix<std::complex<double>> finalMatrix = newMatrix * (*complexRateEigen.cDiagRightMatrix) ;

  // some floating precision error so use EXPECT_NEAR, also we only compare the real counterparts because
  // the imaginary should be close to 0 (very close)
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      EXPECT_NEAR(finalMatrix(i,j).real(), testMatrix(i,j), 0.0001);
    }
  }
}