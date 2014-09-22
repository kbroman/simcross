IntegerMatrix cpp_convert2geno(List xodat, NumericVector map);
IntegerVector convertchr2geno(List chr, NumericVector map);
IntegerMatrix combine_mat_and_pat_geno(IntegerMatrix matmatrix, IntegerMatrix patmatrix, int max_geno);
IntegerMatrix combine_mat_and_pat_geno_wfounders(IntegerMatrix matmatrix, IntegerMatrix patmatrix,
                                                 IntegerMatrix founder_geno);
