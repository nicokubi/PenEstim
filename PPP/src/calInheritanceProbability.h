//
// Created by Gavin on 02-Apr-20.
//

#ifndef CALINHERITANCEPROBABILITY_H
#define CALINHERITANCEPROBABILITY_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Inheritance probability from genotypes
//' 
//' Returns the inheritance probability tensor given the matrix of possible genotypes
//' 
//' @details
//' The base cube `M` defined in the body of `calInheritanceProb` is used to 
//' calculate the inheritance probabilities for each possible genotype. 
//' 
//' The first dimension corresponds to the father, the second dimension to the 
//' mother, and the third dimension to the child, such that `M` should be 
//' indexed as `M(father_index, mother_index, child_index)`. 
//' 
//' Along each dimension, the first index (`0`) corresponds to 
//' wildtype/noncarrier, i.e. neither allele has the mutation at this locus; 
//' the second index (`1`) corresponds to heterozygous carrier, i.e. one of 
//' the alleles is mutated and the other is not; and the third index (`2`) 
//' corresponds to homozygous carrier, i.e. both alleles carry the mutation. 
//' 
//' Each element corresponds to the probability that the child has the 
//' genotype represented by the third dimension, given that their parents have 
//' the genotypes represented by the first and second dimensions. 
//' 
//' Examples: 
//' 
//' `M(1, 0, 2) == 0`: Father is heterozygous carrier, mother is wildtype, and 
//' child is homozygous carrier. This is an impossible configuration, since 
//' there is no way for the child to inherit a mutated allele from their 
//' mother, so the probability is `0`. 
//' 
//' `M(0, 1, 0) == 0.5`: Father is wildtype, mother is heterozygous carrier, 
//' and child is wildtype. The child always inherits a wildtype allele from 
//' their father and has a 50/50 chance of inheriting a wildtype allele from 
//' their mother, so there is a `0.5` probability of being wildtype. 
//' 
//' `M(0, 2, 1) == 1`: Father is wildtype, mother is homozygous carrier, and 
//' child is heterozygous carrier. The child always inherits a wildtype 
//' allele from their father and always inherits a mutated allele from their 
//' mother, so the child must be heterozygous, hence a probability of `1`. 
//' 
//' @param possGenotypes matrix called possGenotypes
//' @return 3D Tensor (cube) giving the transition/inheritance probabilities
// [[Rcpp::export(name=".calInheritanceProb")]]
arma::Cube<double> calInheritanceProb(const arma::Mat<int>& possGenotypes);

#endif //CALINHERITANCEPROBABILITY_H
