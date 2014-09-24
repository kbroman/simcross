
all: vignettes
.PHONY: vignettes

VIGNETTES = assets/vignettes/simcross.html
vignettes: ${VIGNETTES}

assets/vignettes/%.html: ../simcross/vignettes/%.Rmd
	R -e 'library(knitr);knit2html("$<", "$@")'

