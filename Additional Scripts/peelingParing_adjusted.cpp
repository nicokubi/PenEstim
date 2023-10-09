//
// Created by Gavin on 02-Apr-20.
// Adjusted to output matrix with the Likelihood of the pedigre by Nicolas in Aug - 23
//

#include "peelingParingHelpers.h"
#include "peelingParing.h"

arma::Row<double> calAntProb(int idx, const arma::Mat<int>& idMatrix, const arma::Mat<int>& mates,
                             const arma::Cube<int>& children, const arma::Mat<double>& prevalence,
                             const arma::Mat<double>& likelihood, const arma::Cube<double>& transition,
                             arma::Mat<double>& antProb,
                             arma::Cube<double>& postProb) {
    
    // Check if we have already calculated the value and return if so
    if (any(antProb.row(idx) > 0)) {
        
        return antProb.row(idx);
    }
    
    // Get motherIdx and fatherIdx
    int motherIdx = getMotherIdx(idx, idMatrix);
    int fatherIdx = getFatherIdx(idx, idMatrix);
    
    // Get the number of possible genotypes considered
    int nGenotypes = prevalence.n_cols;

    // Initialise result matrix
    // This represents the two dimensions of mother and father
    // The third dimension will represent the children
    // Reference: Fernando et al 1993: Equation (3) line 4
    arma::Mat<double> childMatrix(nGenotypes, nGenotypes, arma::fill::ones);

    // We begin with the inner most calculation and work outwards

    // a. Get siblings (shared parents), ref: product over j in C_{m,f}, j =/= i
    arma::Col<int> sibsIdxs = getSibsIdxs(idx, idMatrix, children);
    
    // Loop over the siblings
    for (arma::uword j = 0; j < sibsIdxs.n_elem; j++) {

        // Get the mates of this particular sibling, ref: product over k in S_j,
        arma::Col<int> matesOfSib = getMatesIdxs(sibsIdxs(j), mates);

        // Intialise vector for the product over k in S_j, j are siblings
        // Copy the initial value as the likelihood of the sibling, ref: g(y_j | u_j)
        arma::Row<double> postProbProd(likelihood.row(sibsIdxs(j)));

        // Perform the product over k in S_j
        for (arma::uword k = 0; k < matesOfSib.n_elem; k++) {

            postProbProd %= calPostProb(sibsIdxs(j), matesOfSib(k), idMatrix, mates, children, prevalence,
                                        likelihood, transition, antProb, postProb);
        }

        // Place holder for multiplication with transition(u_j | u_m, u_f)
        // Copy transition as the initial value, ref: Equation (3) line 4
        arma::Cube<double> childCub(transition);

        // Check dimensions match
        if (postProbProd.n_elem != childCub.n_slices) {
            Rcpp::stop("Dimension mismatch. Couldn't multiply cube with N slices with vector of N elements...");
        }

        for (arma::uword i = 0; i < postProbProd.n_elem; i++) {
            // Multiply each slice by scalar postProbProd(i)
            // That is, each element of the matrix S.slice(i) is multiplied by postProbProd(i)
            // Ref: Equation (3) line 4: tr(u_j | u_m, u_f) * g(y_j | u_j) * prod_{k in S_j} p_{jk} (u_j)
            childCub.slice(i) *= postProbProd(i);
        }

        // Sum childCub over the slices and multiply it to the childMatrix
        // Ref: Equation (3) line 4, %= is the product over j in C_{mf}, j =/= i, and sum is over u_j
        childMatrix %= sum(childCub, 2);
    }

    // b. Mother and father loops
    // Deal with the case when there is only one parent as well
    // Case when there is a mother but no father
    
    // Now take care of mother's loop, ref: Equation (3) line 1 (line 1 and 2 interchangeable)
    // First multiply by anterior * likelihood, ref: a_m(u_m) * g(y_m | u_m)
    if (motherIdx != -999) {
        // There is a mother but no father
        childMatrix.each_col() %= (likelihood.row(motherIdx) % 
                                  calAntProb(motherIdx, idMatrix, mates, children, prevalence, likelihood, transition,
                                            antProb, postProb)).t();
        
    } else if (fatherIdx != -999) {
        // Use the idx's prevalences for the missing parent
        // There is no mother but there is a father
        childMatrix.each_col() %= (prevalence.row(idx)).t();
    }
    
    // c. Get mates of mother who are not the father
    arma::Col<int> mothersOtherMatesIdxs = getOtherMatesIdxs(motherIdx, fatherIdx, mates);
    
    // If there are other mates of the mother, multiply each row
    // Ref: multiply childMatrix by product of posteriors over mother's mates who are not the father
    // Ref: product over j in S_m, j =/= f of p_{mj} (u_m)
    for (arma::uword m = 0; m < mothersOtherMatesIdxs.n_elem; m++) {

        childMatrix.each_col() %= calPostProb(motherIdx, mothersOtherMatesIdxs(m), idMatrix, mates, children,
                                              prevalence, likelihood, transition, antProb, postProb).t();
    }
    
    // Deal with the case when there is only one parent as well
    // Case when there is a father but no mother
    
    // Now take care of father's loop, ref: Equation (3) line 2 (line 1 and 2 interchangeable)
    // First multiply by anterior * likelihood, ref: a_f(u_f) * g(y_f | u_f)
    if (fatherIdx != -999) {
        // There is a father
        childMatrix.each_row() %= likelihood.row(fatherIdx) %
                                  calAntProb(fatherIdx, idMatrix, mates, children, prevalence, likelihood, transition,
                                               antProb, postProb);
        
    } else if (motherIdx != -999) {
        // Use the idx's prevalences for the missing parent
        // There is no father but there is a mother
        childMatrix.each_row() %= prevalence.row(idx);
    }
    
    // d. Get mates of father who are not the mother
    arma::Col<int> fathersOtherMatesIdxs = getOtherMatesIdxs(fatherIdx, motherIdx, mates);

    // If there are other mates of the father, multiply each column
    // Ref: multiply childMatrix by product of posteriors over father's mates who are not the mother
    // Ref: product over j in S_f, j =/= m of p_{fj} (u_f)
    for (arma::uword f = 0; f < fathersOtherMatesIdxs.n_elem; f++) {

        childMatrix.each_row() %= calPostProb(fatherIdx, fathersOtherMatesIdxs(f), idMatrix, mates, children,
                                              prevalence, likelihood, transition, antProb, postProb);
    }

    // Initialise innerCube with the transition matrix by copy
    // Ref: Equation (3) line 3, tr(u_i | u_m, u_f)
    arma::Cube<double> innerCube(transition);

    // Multiply each slice by the childMatrix, ref: x sign in Equation (3) line 4
    innerCube.each_slice() %= childMatrix;

    // Sum the innerCube over the mother and father to get row vector (implicit conversion by copy)
    // Ref: sum over u_m and u_f in Equation (3) line 1 and 2
    arma::Row<double> antProbRow = sum(sum(innerCube, 0), 1);

    // Set by reference for later
    antProb.row(idx) = antProbRow;

    return antProbRow;
}

arma::Row<double>
calPostProb(int idx, int mateIdx, const arma::Mat<int>& idMatrix, const arma::Mat<int>& mates,
            const arma::Cube<int>& children, const arma::Mat<double>& prevalence,
            const arma::Mat<double>& likelihood, const arma::Cube<double>& transition,
            arma::Mat<double>& antProb, arma::Cube<double>& postProb) {

    // Check that we aren't calculating this for someone who doesn't exist in the pedigree
    if (idx == -999 || mateIdx == -999) {
        Rcpp::stop("idx = -999 or mateIdx = -999 ... Missing idx or mateIdx");
    }
    
    // Check if we have already calculated the value and return if so
    arma::Row<double> tube = postProb(arma::span(idx), arma::span(mateIdx), arma::span::all);

    if (any(tube) > 0) {
        
        return tube;
    }

    // Get the number of possible genotypes considered
    int nGenotypes = prevalence.n_cols;

    // Initialise result matrix
    // This represents the two dimensions of mother and father
    // The third dimension will represent the children
    // Reference: Equation (4) line 2, product over k in C_{ij}
    arma::Mat<double> childMatrix(nGenotypes, nGenotypes, arma::fill::ones);

    // a. Get children of idx and mateIdx
    arma::Col<int> childrenIdxs = getChildrenIdxs(idx, mateIdx, children);

    // Loop over children, ref: Equation (4) line 2, product over k in C_{ij}
    for (arma::uword k = 0; k < childrenIdxs.n_elem; k++) {
        
        // Get the mates of this particular child, ref: product over l in S_k
        // Note that if there are no real mates, getMatesIdxs will give the idx of the input
        arma::Col<int> matesOfChild = getMatesIdxs(childrenIdxs(k), mates);

        // Initialise vector for postProbProduct result
        // Ref: Equation (4) line 2, prod_{l in mates of k} p_{kl} (u_k), where k are children
        // Initialise this with the likelihood, which is g(y_k | u_k)
        arma::Row<double> postProbProd(likelihood.row(childrenIdxs(k)));

        // Ref: Equation (4) line 2, calculate product over l for l in S_k
        for (arma::uword l = 0; l < matesOfChild.n_elem; l++) {
            
            postProbProd %= calPostProb(childrenIdxs(k), matesOfChild(l), idMatrix, mates, children, prevalence,
                                        likelihood, transition, antProb, postProb);

        }

        // Place holder for multiplication with transition(u_k | u_i, u_j)
        // Give it a copy of transition, ref: Equation (4) line (2)
        arma::Cube<double> childCub(transition);

        // Check dimensions match
        if (postProbProd.n_elem != childCub.n_slices) {
            Rcpp::stop("Dimension mismatch. Couldn't multiply cube with N slices with vector of N elements...");
        }

        for (arma::uword i = 0; i < postProbProd.n_elem; i++) {
            // Multiply each slice by scalar postProbProd(i)
            // That is, each element of the matrix S.slice(i) is multiplied by postProbProd(i)
            childCub.slice(i) *= postProbProd(i);
        }

        // Sum childCub over the slices and multiply childMatrix elementwise
        // Ref: Equation (4) line 2, sum over u_k inside and product over k in C_{ij} outside
        childMatrix %= sum(childCub, 2);

    }

    // Initialise mate matrix as childMatrix to continue calculation
    // Ref: Equation (4) line 1, use childMatrix, which represents the whole of line 2
    arma::Mat<double> mateMatrix(childMatrix);

    // Deal with the case when there is a single parent
    // If idx and mateIdx are different, with childrenIdxs.n_elem > 0
    // then initialise the mate row
    // Else use the prevalence of the index for the missing parent
    
    if (idx != mateIdx && childrenIdxs.n_elem > 0) {
        // Both parents are in the pedigree and they have children
        // Initialise row mateRow to multiply the otherMatesIdxs matrix
        // Ref: Equation (4) line 1, a_j(u_j) * g(y_j | u_j)
        arma::Row<double> mateRow = likelihood.row(mateIdx) %
            calAntProb(mateIdx, idMatrix, mates, children, prevalence, likelihood, transition,
                       antProb, postProb);
        
        // Sweep mate matrix by mateRow, ref: Equation (4) line 1, a_j(u_j) * g(y_j | u_j) * prod over k in C_{ij}
        mateMatrix.each_col() %= mateRow.t();
        
    } else if (idx == mateIdx && childrenIdxs.n_elem > 0) {
        // Sweep mate matrix by the prevalence in the case that there is a single parent
        mateMatrix.each_col() %= prevalence.row(childrenIdxs(0)).t();
    }

    // b. Get mateIdxs mates other than idx
    // Ref: get the product over k in S_j, k =/= i, p_{jk} (u_j)
    arma::Col<int> otherMatesIdxs = getOtherMatesIdxs(mateIdx, idx, mates);

    // Ref: calculate product over k in S_j, k =/= i, p_{jk} (u_j)
    for (arma::uword k = 0; k < otherMatesIdxs.n_elem; k++) {

        mateMatrix.each_row() %= calPostProb(mateIdx, otherMatesIdxs(k), idMatrix, mates, children, prevalence,
                                             likelihood, transition, antProb, postProb);

    }

    // Sum the mateMatrix over the mates (columns) row vector, ref: sum over u_j in Equation (4) line 1
    // Set by reference for later
    postProb.tube(idx, mateIdx) = sum(mateMatrix, 0);

    return postProb.tube(idx, mateIdx);
}

arma::Mat<double>
peelingParing(arma::Col<int> probandIdxs, const arma::Mat<int>& idMatrix, const arma::Mat<int>& mates,
              const arma::Cube<int>& children,
              const arma::Mat<double>& prevalence, const arma::Mat<double>& likelihood,
              const arma::Cube<double>& transition, arma::Mat<double>& antProb, arma::Cube<double>& postProb) {

    // Get the number of possible genotypes considered
    int nGenotypes = prevalence.n_cols;

    // Initialise the table to be returned
    arma::Mat<double> posteriorProbs(probandIdxs.n_elem, nGenotypes, arma::fill::zeros);

    // Initialise row of ones outside the loop to avoid repeated memory allocation
    arma::Row<double> postProbProd(nGenotypes, arma::fill::ones);
    
    for (arma::uword i = 0; i < probandIdxs.n_elem; i++) {

        // Reset postProbProd to ones at the beginning of each iteration
        postProbProd.fill(1.0);

        // Get the mates of probandIdx and calculate the product of the postProb
        arma::Col<int> matesOfProbandIdxs = getMatesIdxs(probandIdxs(i), mates);

        for (arma::uword j = 0; j < matesOfProbandIdxs.n_elem; j++) {
            postProbProd %= calPostProb(probandIdxs(i), matesOfProbandIdxs(j), idMatrix, mates, children,
                                        prevalence, likelihood, transition, antProb, postProb);
        }

        // Ref: Equation (1)
        posteriorProbs.row(i) = calAntProb(probandIdxs(i), idMatrix, mates, children, prevalence, likelihood,
                                           transition, antProb, postProb)
                                % likelihood.row(probandIdxs(i))
                                % postProbProd;
    }

    // Calculate the normalization constants (row sums) as a column vector
    arma::Col<double> normalizationConstants = sum(posteriorProbs, 1);
    
    // Construct the output matrix by replicating the column vector across columns
    arma::Mat<double> outputMatrix = repmat(normalizationConstants, 1, nGenotypes);
    
    return outputMatrix;
}
