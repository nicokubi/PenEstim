//
// Created by Gavin on 14-Apr-20.
//

#ifndef CALPREVALENCES_H
#define CALPREVALENCES_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute prevalences for each admissible Genotype
//' 
//' Returns the inheritance probability tensor given the matrix of possible genotypes
//'
//' @param margPrevs the marginal prevalence matrix of each mutation in each pedigree member
//' @param possGenotypes matrix called possGenotypes
//' @return Matrix giving the prevalences for each pedigree member
// [[Rcpp::export(name = ".calPrevalences")]]
arma::Mat<double> calPrevalences(const arma::Mat<double>& margPrevs, const arma::Mat<int>& possGenotypes);

#endif //CALPREVALENCES_H
