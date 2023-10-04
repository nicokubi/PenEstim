//
// Created by Gavin on 02-Apr-20.
//

#include "calInheritanceProbability.h"

arma::Cube<double> calInheritanceProb(const arma::Mat<int>& possGenotypes) {
    arma::uword nGenotypes = possGenotypes.n_rows;
    arma::uword nGenes = possGenotypes.n_cols;

    // Initialise the cube to be returned
    arma::Cube<double> transition(nGenotypes, nGenotypes, nGenotypes, arma::fill::zeros);

    // Initialise 3x3x3 base cube
    // The first dimension corresponds to the father, the second dimension to 
    // the mother, and the third dimension to the child, such that M should be 
    // indexed as M(father_index, mother_index, child_index). 
    // Along each dimension, the first index (0) corresponds to 
    // wildtype/noncarrier, i.e. neither allele has the mutation at this 
    // locus; the second index (1) corresponds to heterozygous carrier, i.e. 
    // one of the alleles is mutated and the other is not; and the third index 
    // (2) corresponds to homozygous carrier, i.e. both alleles carry the 
    // mutation. 
    // Each element corresponds to the probability that the child has the 
    // genotype represented by the third dimension, given that their parents 
    // have the genotypes represented by the first and second dimensions. 
    // Examples: 
    // M(1, 0, 2) == 0: Father is heterozygous carrier, mother is wildtype, 
    // and child is homozygous carrier. This is an impossible configuration, 
    // since there is no way for the child to inherit a mutated allele from 
    // their mother, so the probability is 0. 
    // M(0, 1, 0) == 0.5: Father is wildtype, mother is heterozygous carrier, 
    // and child is wildtype. The child always inherits a wildtype allele from 
    // their father and has a 50/50 chance of inheriting a wildtype allele 
    // from their mother, so there is a 0.5 probability of being wildtype. 
    // M(0, 2, 1) == 1: Father is wildtype, mother is homozygous carrier, and 
    // child is heterozygous carrier. The child always inherits a wildtype 
    // allele from their father and always inherits a mutated allele from 
    // their mother, so the child must be heterozygous, hence a probability of 
    // 1. 
    arma::Cube<double> M(3, 3, 3, arma::fill::zeros);
    M.slice(0) = {{1,   0.5,  0},
                  {0.5, 0.25, 0},
                  {0,   0,    0}};
    M.slice(1) = {{0,   0.5, 1},
                  {0.5, 0.5, 0.5},
                  {1,   0.5, 0}};
    M.slice(2) = {{0, 0,    0},
                  {0, 0.25, 0.5},
                  {0, 0.5,  1}};

    for (arma::uword i = 0; i < nGenotypes; i++) {

        arma::Row<int> u_f = possGenotypes.row(i);  // Rows of possGenotypes represent genotypes

        for (arma::uword j = 0; j <= i; j++) {

            arma::Row<int> u_m = possGenotypes.row(j);

            for (arma::uword l = 0; l < nGenotypes; l++) {

                arma::Row<int> u = possGenotypes.row(l);

                double prod = 1; // Placeholder

                for (arma::uword k = 0; k < nGenes; k++) {

                    double prob = M(u_f(k), u_m(k), u(k));

                    if (prob == 0) {
                        prod = 0;
                        break;
                    } else {
                        prod = prod * prob;
                    }
                }
                transition(i, j, l) = prod;
            }
            transition.tube(j, i) = transition.tube(i, j);
        }
    }

    return transition;
}
