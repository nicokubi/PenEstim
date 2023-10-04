//
// Created by Gavin on 02-Apr-20.
//

#ifndef PEELINGPARINGHELPERS_H
#define PEELINGPARINGHELPERS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' getFatherIdx
//' 
//' Gets the index number of the Father of proband
//' 
//' Columns in idMatrix are 0, 1, 2 (ID, MotherID, FatherID)
//' @param idx The idx number of the proband in the ID column (C++ row number)
//' @param idMatrix The matrix representing the ID numbers of the ID, MotherID and FatherID
//' @return Father's index (row number)
// [[Rcpp::export(name=".getFatherIdx")]]
int getFatherIdx(int idx, const arma::Mat<int>& idMatrix);

//' getMotherIdx
//' 
//' Gets the index number of the mother of proband
//' 
//' Columns in idMatrix are 0, 1, 2 (ID, MotherID, FatherID)
//' @param idx The idx number of the proband in the ID column (C++ row number)
//' @param idMatrix The matrix representing the ID numbers of the ID, MotherID and FatherID
//' @return Mother's index (row number)
// [[Rcpp::export(".getMotherIdx")]]
int getMotherIdx(int idx, const arma::Mat<int>& idMatrix);

//' getSibsIdxs
//' 
//' Gets the index numbers of possible siblings of proband
//' 
//' Columns in idMatrix are 0, 1, 2 (ID, MotherID, FatherID)
//' @param idx The idx number of the proband in the ID column (C++ row number)
//' @param idMatrix The matrix representing the ID numbers of the ID, MotherID and FatherID
//' @param children The `cube' specifying the children
//' @return Indexes of siblings
// [[Rcpp::export(".getSibsIdxs")]]
arma::Col<int> getSibsIdxs(int idx, const arma::Mat<int>& idMatrix, const arma::Cube<int>& children);

//' getMatesIdxs
//' 
//' Gets the index numbers of possible mates of proband
//' @param idx The idx number of the proband (C++ row number)
//' @param mates The `matrix' specifying the mates
//' @return Idxs of mates
// [[Rcpp::export(".getMatesIdxs")]]
arma::Col<int> getMatesIdxs(int idx, const arma::Mat<int>& mates);

//' getOtherMatesIdxs
//'
//' Gets the index numbers of possible mates of idx who are not mateIdx
//' @param idx The idx number of the proband (C++ row number)
//' @param mateIdx The idx number of the mate of the proband
//' @param mates The `matrix' specifying the mates
//' @return Idxs of other mates
// [[Rcpp::export(".getOtherMatesIdxs")]]
arma::Col<int> getOtherMatesIdxs(int idx, int mateIdx, const arma::Mat<int>& mates);

//' getChildrenIdxs
//' 
//' Gets the index numbers of possible children between idx and mateIDX
//' @param idx The idx number of the proband (C++ row number)
//' @param mateIdx The other parent
//' @param children The `cube' specifying the children
//' @return Idxs of chidlren
// [[Rcpp::export(".getChildrenIdxs")]]
arma::Col<int> getChildrenIdxs(int idx, int mateIdx, const arma::Cube<int>& children);

#endif //PEELINGPARINGHELPERS_H
