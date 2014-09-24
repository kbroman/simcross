---
layout: page
title: R/simcross
tagline: simulation of general experimental crosses
---

[R/simcross](https://github.com/kbroman/simcross) is an
[R](http://www.r-project.org) package to simulate general experimental
crosses, including advanced intercross lines and the diversity outcross.
It is intended to be a companion to [R/qtl](http://www.rqtl.org).

### Installation

R/simcross is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/simcross from its
[GitHub repository](http://github.com/kbroman/simcross). You first need to
install the [devtools](https://github.com/hadley/devtools) package.

    install.packages("devtools")

Then install R/simcross using the `install_github` function in the
[devtools](http://github.com/hadley/devtools) package.

    library(devtools)
    install_github("kbroman/simcross")

### Vignette

A vignette describing the use of the package is available from within
R (and [also here](assets/vignettes/simcross.html). Load the package
and then use the `vignette` function.

    library(simcross)
    vignette("simcross", package="simcross")

### License

The R/simcross package is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>

---

Sources on [github](http://github.com):

- The [source for the package](https://github.com/kbroman/simcross/tree/master)
- The [source for the website](https://github.com/kbroman/simcross/tree/gh-pages)
