## simcross 0.2-29 (2018-08-16)

### Minor changes

- Added dataset, mouseL_cox and mouseL_mgi, with mouse chromosome
  lengths in cM.

- In C++, avoid calling runif() with n=0, and require Rcpp 0.12.17
  (regarding handling of random number generation)

- Byte-compile the R code


## simcross 0.2-25 (2017-11-21)

### New features

- Added `palette` argument to `CCcolors()`. Default is to return my
  revised CC colors. `CCcolors("orig")` gives the original, official
  colors.

- The `chrlength` argument for `plot_ind()` can be a vector of length
  2, in which case the two chromosomes will be different lengths but
  aligned at the top. This is to allow illustratation of the X and Y
  chromosomes in males.


## simcross 0.2-22 (2017-04-05)

### New features

- Revised sim_ril_pedigree to work for 2^n-way RIL for any n>0.
