//
// Created by Gavin on 14-Apr-20.
//

#include "calPrevalences.h"

arma::Mat<double> calPrevalences(const arma::Mat<double>& margPrevs, const arma::Mat<int>& possGenotypes) {
  // Initialise prevalences
  arma::Mat<double> prevalences(margPrevs.n_rows, possGenotypes.n_rows, arma::fill::zeros);
  
  // Number of family members in pedigree
  arma::uword N = margPrevs.n_rows;
  
  // Number of mutations
  arma::uword K = margPrevs.n_cols;
  
  // Initialise M matrix
  arma::Cube<double> M(N, K, 3, arma::fill::zeros);
  
  // Initialise `ones' matrix of the same size as margPrevs
  arma::Mat<double> Ones(margPrevs.n_rows, margPrevs.n_cols, arma::fill::ones);
  
  // Implement formula for M: M = (1-margPrevs)^2, 2 * margPrevs * (1-margPrevs), margPrevs^2 elementwise
  M.slice(0) = pow(Ones - margPrevs, 2);
  M.slice(1) = 2 * margPrevs % (Ones - margPrevs);
  M.slice(2) = pow(margPrevs, 2);
  
  // Loop over number of family in pedigree
  for (arma::uword n = 0; n < N; n++) {
    // Loop over possible genotypes
    for (arma::uword i = 0; i < possGenotypes.n_rows; i++) {
      // Start product
      double prev = 1;
      
      // Multiply over the M slices
      for (arma::uword j = 0; j < K; j++) {
        prev *= M(n, j, possGenotypes(i, j));
      }
      // Set value in the prevalence matrix
      prevalences(n, i) = prev;
    }
  }
  
  return prevalences;
}
