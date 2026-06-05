## [R/simcross](https://kbroman.org/simcross/) <a href="https://kbroman.org/simcross/"><img src="https://kbroman.org/simcross/assets/pics/simcross_logo.png" align="right" height="138" alt="R/simcross logo"/></a>

[![R-CMD-check](https://github.com/kbroman/simcross/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kbroman/simcross/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/simcross)](https://cran.r-project.org/package=simcross)
[![r-universe badge](https://kbroman.r-universe.dev/simcross/badges/version)](https://kbroman.r-universe.dev/simcross)
[![zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4032914.svg)](https://doi.org/10.5281/zenodo.4032914)

[Karl W Broman](https://kbroman.org)

---

R/simcross is an [R](https://www.r-project.org) package for simulating
and plotting general experimental crosses.

---

### Installation

You can install R/simcross from [CRAN](https://cran.r-project.org):

```r
install.packages("simcross")
```

Alternatively, install it from [R
universe](https://kbroman.r-universe.dev):

```r
install.packages("simcross", repos=c("https://kbroman.r-universe.dev",
                                     "https://cloud.r-project.org"))
```

Or use [remotes](https://remotes.r-lib.org) to install it from its GitHub source:

```r
install.packages("remotes")
remotes::install_github("kbroman/simcross")
```

---

### Vignette

A vignette describing the use of the package is available from within
R (and [also here](https://kbroman.org/simcross/assets/vignettes/simcross.html)). Load the package
and then use the `vignette` function.

```r
library(simcross)
vignette("simcross", package="simcross")
```

---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>
