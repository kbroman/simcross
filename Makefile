all: vignettes
.PHONY: vignettes

VIGNETTES = assets/vignettes/simcross.html
vignettes: ${VIGNETTES}

assets/vignettes/%.html: ../simcross/vignettes/%.Rmd
	cd $(<D); \
	R -e "rmarkdown::render('$(<F)')"; \
	mv $(@F) ../../Web/$(@D)/
