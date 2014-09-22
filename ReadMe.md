## R/simcross

[![Build Status](https://travis-ci.org/kbroman/simcross.png?branch=master)](https://travis-ci.org/kbroman/simcross)

copyright (c) 2014, [Karl W Broman](http://www.biostat.wisc.edu/~kbroman)

---

R/simcross is an [R](http://www.r-project.org) package for simulating
and plotting general experimental crosses.

---

### Installation

R/simcross is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/simcross from its
[GitHub repository](http://github.com/kbroman/simcross). You first need to
install the [devtools](https://github.com/hadley/devtools) package.

```S
install.packages("devtools")
```

Then install R/simcross using the `install_github` function in the
[devtools](http://github.com/hadley/devtools) package.

```S
library(devtools)
install_github("kbroman/simcross")
```

---

### Vignette

A vignette describing the use of the package is available from within
R. Load the package and then use the `vignette` function.

```S
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
<http://www.r-project.org/Licenses/GPL-3>
