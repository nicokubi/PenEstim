//
// Created by Gavin on 02-Apr-20.
//

#ifndef PEELINGPARING_H
#define PEELINGPARING_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate Anterior Probability
//' 
//' Calculates the anterior probability of idx
//' 
//' We assume that we have already checked for loops and that there are no loops
//' 
//' @param idx The idx number of the proband (C++ row number)
//' @param idMatrix The matrix representing the idx numbers of the ID, MotherID and FatherID
//' @param mates The `matrix' specifying the mates
//' @param children The `cube' specifying the children
//' @param prevalence The `matrix' specifying the prevalences
//' @param likelihood The `matrix' specifying the likelihoods (penetrances)
//' @param transition The `cube' specifying the transition probabilities from mother and father to child
//' @param antProb A placeholder to memoise the recursive function output of calAntProb, size #Members x #Genotypes
//' @param postProb A placeholder to memoise the recursive function output of calPostProb, size #Members x #Members x #Genotypes
//' @return Values of the anterior probability of idx having possible genotypes
// [[Rcpp::export(name = ".calAntProb")]]
arma::Row<double> calAntProb(int idx, const arma::Mat<int>& idMatrix, const arma::Mat<int>& mates,
                             const arma::Cube<int>& children, const arma::Mat<double>& prevalence,
                             const arma::Mat<double>& likelihood, const arma::Cube<double>& transition,
                             arma::Mat<double>& antProb,
                             arma::Cube<double>& postProb);

//' Calculate Posterior Probability
//' 
//' Calculates the posterior probability of idx and mateIdx
//' 
//' We assume that we have already checked for loops and that there are no loops
//' @param idx The idx number of the proband (C++ row number)
//' @param mateIdx The idx number of the mate of `idx' (C++ row number)
//' @param idMatrix The matrix representing the idx numbers of the ID, MotherID and FatherID
//' @param mates The `matrix' specifying the mates
//' @param children The `cube' specifying the children
//' @param prevalence The `matrix' specifying the prevalences
//' @param likelihood The `matrix' specifying the likelihoods (penetrances)
//' @param transition The `cube' specifying the transition probabilities from mother and father to child
//' @param antProb A placeholder to memoise the recursive function output of calAntProb, size #Members x #Genotypes
//' @param postProb A placeholder to memoise the recursive function output of calPostProb, size #Members x #Members x #Genotypes
//' @return Values of the posterior probability of idx and mateIdx having possible genotypes
// [[Rcpp::export(name = ".calPostProb")]]
arma::Row<double>
calPostProb(int idx, int mateIdx, const arma::Mat<int>& idMatrix, const arma::Mat<int>& mates,
            const arma::Cube<int>& children, const arma::Mat<double>& prevalence,
            const arma::Mat<double>& likelihood, const arma::Cube<double>& transition,
            arma::Mat<double>& antProb, arma::Cube<double>& postProb);

//' Peeling and Paring
//' 
//' Calculates the posterior probabilities of genotypes of probands given inputs
//' 
//' Calls calAntProb and calPostProb recursively
//' NB: Require that NAs in antProb and postProb are set to 0
//' @param probandIdxs The idx numbers of the probands in the ID column
//' @param idMatrix The matrix representing the idx numbers of the ID, MotherID and FatherID
//' @param mates The `matrix' specifying the mates
//' @param children The `cube' specifying the children
//' @param prevalence The `matrix' specifying the prevalences
//' @param likelihood The `matrix' specifying the likelihoods (penetrances)
//' @param transition The `cube' specifying the transition probabilities from mother and father to child
//' @param antProb A placeholder to memoise the recursive function output of calAntProb
//' @param postProb A placeholder to memoise the recursive function output of calPostProb
//' @return Matrix of the posterior probabilities of ID having certain genotypes
// [[Rcpp::export(".peelingParing")]]
arma::Mat<double>
peelingParing(arma::Col<int> probandIdxs, const arma::Mat<int>& idMatrix, const arma::Mat<int>& mates,
              const arma::Cube<int>& children,
              const arma::Mat<double>& prevalence, const arma::Mat<double>& likelihood,
              const arma::Cube<double>& transition, arma::Mat<double>& antProb, arma::Cube<double>& postProb);

#endif //PEELINGPARING_H
