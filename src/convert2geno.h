IntegerMatrix fromR_convert2geno(const List xodat, const NumericVector map);
IntegerVector convertchr2geno(const List chr, const NumericVector map);
IntegerMatrix combine_mat_and_pat_geno(const IntegerMatrix matmatrix, const IntegerMatrix patmatrix, const int max_geno);
IntegerMatrix combine_mat_and_pat_geno_wfounders(const IntegerMatrix matmatrix, const IntegerMatrix patmatrix,
                                                 const IntegerMatrix founder_geno);
IntegerVector fromR_convert2genoarray(const List xodat, const NumericVector map);
