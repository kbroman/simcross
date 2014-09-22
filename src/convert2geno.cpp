#include <Rcpp.h>
using namespace Rcpp;

#include "convert2geno.h"

// [[Rcpp::export]]
IntegerMatrix cpp_convert2geno(List xodat, NumericVector map)
{
    int n_ind = xodat.size();
    int n_mar = map.size();
    IntegerMatrix matmatrix(n_mar, n_ind), patmatrix(n_mar, n_ind);

    for(int i=0; i<n_ind; i++) {
        List ind = xodat[i];
        List mat = ind[0];
        List pat = ind[1];

        // convert mat and pat chr to genotype matrices
        IntegerVector matdata = convertchr2geno(mat, map);
        std::copy(matdata.begin(), matdata.end(), matmatrix.begin()+i*n_mar);
        IntegerVector patdata = convertchr2geno(pat, map);
        std::copy(patdata.begin(), patdata.end(), patmatrix.begin()+i*n_mar);
    }

    // find maximum genotype
    int max_geno_mat = max(matmatrix);
    int max_geno_pat = max(patmatrix);
    int max_geno = max_geno_mat > max_geno_pat ? max_geno_mat : max_geno_pat;

    IntegerVector result_vector = combine_mat_and_pat_geno(matmatrix, patmatrix, max_geno);
    IntegerMatrix result(n_mar, n_ind);
    std::copy(result_vector.begin(), result_vector.end(), result.begin());

    return result;
}


IntegerVector convertchr2geno(List chr, NumericVector map)
{
    IntegerVector alleles = chr[0];
    NumericVector locations = chr[1];

    int n_mar = map.size();
    int n_loc = locations.size();

    IntegerVector output(n_mar);

    // at each marker, find closest position to right and take allele there
    for(int i=0; i<n_mar; i++) {
        for(int j=0; j<n_loc; j++) {
            if(locations[j] >= map[i]) {
                output[i] = alleles[j];
                break;
            }
        }
    }

    return output;
}

IntegerVector combine_mat_and_pat_geno(IntegerVector matmatrix, IntegerVector patmatrix, int max_geno)
{
    if(max_geno == 2) {
        for(int i=0; i<matmatrix.size(); i++)
            matmatrix[i] = matmatrix[i] + patmatrix[i] - 1;
    }
    else {  // multi-allele case
        for(int i=0; i<matmatrix.size(); i++)
            matmatrix[i] = (1 << (matmatrix[i]-1)) + (1 << (patmatrix[i] - 1));
    }
    
    return matmatrix;
}
