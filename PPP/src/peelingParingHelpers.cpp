//
// Created by Gavin on 02-Apr-20.
//

#include "peelingParingHelpers.h"

int getFatherIdx(int idx, const arma::Mat<int>& idMatrix) {
    
    // Check if person is founder
    if (idMatrix(idx, 2) == -999) { return -999; } else {

    // Find index
    arma::uvec ID_bool = find(idMatrix.col(0) == idMatrix(idx, 2));

    return ID_bool(0);
    
    }
}

int getMotherIdx(int idx, const arma::Mat<int>& idMatrix) {
    
    // Check if person is founder
    if (idMatrix(idx, 1) == -999) { return -999; } else {

    // Find index
    arma::uvec ID_bool = find(idMatrix.col(0) == idMatrix(idx, 1));

    return ID_bool(0);
    
    }
}

arma::Col<int> getSibsIdxs(int idx, const arma::Mat<int>& idMatrix, const arma::Cube<int>& children) {
    
    // Get the mother and father Idxs
    int motherIdx = getMotherIdx(idx, idMatrix);
    int fatherIdx = getFatherIdx(idx, idMatrix);
    
    // Deal with the case when either the mother or father are missing
    // The pp.childMatrix function duplicates the single parent
    // So we just need to retrieve them from the same index if one is missing
    
    // Placeholder for empty vector
    arma::Col<int> sibsIdxs(0, arma::fill::zeros);
    
    // Get the possible siblings through the .tube(row, col) function
    // Initialise vector for the possible siblings
    arma::Col<int> possSibs(idMatrix.n_rows, arma::fill::zeros);
    
    if (motherIdx == -999 && fatherIdx != -999) {
        // Has a father but missing mother
        possSibs = children.tube(fatherIdx, fatherIdx);
        
    } else if (motherIdx != -999 && fatherIdx == -999) {
        // Has a mother but missing father
        possSibs = children.tube(motherIdx, motherIdx);
        
    } else if (motherIdx != -999 && fatherIdx != -999) {
        // Both parents present
        possSibs = children.tube(motherIdx, fatherIdx);
        
    } else if (motherIdx == -999 && fatherIdx == -999) {
        // This person is a founder and so return an empty vector
        return sibsIdxs;
    }
    
    // There will be 1 where there are siblings and 0 otherwise
    arma::uvec possSibsIdxs = find(possSibs == 1);
    
    // If there are no siblings, return empty vector
    if (possSibsIdxs.n_elem == 1 || possSibsIdxs.n_elem == 0) { return sibsIdxs; } else {
        
        // Number of siblings will be 1 less than the possible siblings
        int nSibs = possSibsIdxs.size() - 1;
        // Resize to nSibs and then set all elements to 0
        sibsIdxs.zeros(nSibs);
        
        // Fill in the rest of the siblings excluding proband
        unsigned int j = 0; // Keeps track of the index to set
        for (arma::uword i = 0; i < possSibsIdxs.size(); i++) {
            if ((int) possSibsIdxs(i) != idx) {
                sibsIdxs(j) = possSibsIdxs(i);
                j++;
            }
        }
        
        return sibsIdxs;
    }
}

arma::Col<int> getMatesIdxs(int idx, const arma::Mat<int>& mates) {
    
    // Note that mates is symmetric so we may just take the probandIdx column
    arma::Col<int> possMates = mates.col(idx);
    arma::Col<int> matesIdxs = arma::conv_to< arma::Col<int> >::from(find(possMates == 1));
    
    // For the purposes of the peeling, if there is no mate (i.e. single parent with children)
    // then put themselves as the mate
    if (matesIdxs.n_elem == 0) {
        matesIdxs.resize(1);
        matesIdxs(0) = idx;
    }
    
    return matesIdxs;
}

arma::Col<int> getOtherMatesIdxs(int idx, int mateIdx, const arma::Mat<int>& mates) {
    
    // Placeholder for empty vector
    arma::Col<int> otherMatesIdxs(0, arma::fill::zeros);
    
    // If idx is missing, i.e. -999, they should have no mates
    if (idx == -999) { return otherMatesIdxs; }
    
    // Get the mates of idx who are not mateIdx
    arma::Col<int> matesIdxs = getMatesIdxs(idx, mates);

    // If there is one or no mates, return empty vector
    if (matesIdxs.n_elem == 1 || matesIdxs.n_elem == 0) { return otherMatesIdxs; }

    // Number of other mates should be one less than the mates
    int nOtherMates = matesIdxs.n_elem - 1;
    // Resize to nOtherMates and then set all elements to 0
    otherMatesIdxs.zeros(nOtherMates);

    // Fill in the otherMatesIDs
    unsigned int j = 0; // Keeps track of the index
    for (arma::uword i = 0; i < matesIdxs.n_elem; i++) {
        // Add to list only if idx nor mateIdx is already listed
        if (matesIdxs(i) != mateIdx && matesIdxs(i) != idx) {
            otherMatesIdxs(j) = matesIdxs(i);
            j++;
        }
    }
    return otherMatesIdxs;
}

arma::Col<int> getChildrenIdxs(int idx, int mateIdx, const arma::Cube<int>& children) {
    
    // Use cube to find possible Children
    // We now allow idx and mateIdx to be the same
    // when there is a missing parent, the children cube will duplicate the parent and so idx and mateIdx 
    // will be the same
    arma::Col<int> possChildren = children.tube(idx, mateIdx);
    arma::uvec childrenIdxs = find(possChildren == 1);

    return arma::conv_to< arma::Col<int> >::from(childrenIdxs);
}
