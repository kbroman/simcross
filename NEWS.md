## simcross 0.2-26 (2018-03-08)

### Minor changes

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
