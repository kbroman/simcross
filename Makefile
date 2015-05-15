all: vignettes
.PHONY: vignettes

VIGNETTES = assets/vignettes/simcross.html
vignettes: ${VIGNETTES}

assets/vignettes/%.html: ../simcross/vignettes/%.Rmd
	cd $(<D); \
	R -e "library(rmarkdown);render('$(<F)', output_format='html_document')"; \
	mv $(@F) ../../Web/$(@D)/
