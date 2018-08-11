## [R/simcross](http://kbroman.org/simcross)

[![Build Status](https://travis-ci.org/kbroman/simcross.svg?branch=master)](https://travis-ci.org/kbroman/simcross)

[Karl W Broman](http://kbroman.org)

---

R/simcross is an [R](https://www.r-project.org) package for simulating
and plotting general experimental crosses.

---

### Installation

You can install R/simcross from its
[GitHub repository](https://github.com/kbroman/simcross). You first need to
install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install R/simcross using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package. (With
`build_vignettes=TRUE`, the vignettes will be built and installed.)

```r
library(devtools)
install_github("kbroman/simcross", build_vignettes=TRUE)
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
