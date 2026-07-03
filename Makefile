all: vignettes
.PHONY: vignettes

VIGNETTES = assets/vignettes/simcross.html
vignettes: ${VIGNETTES}

assets/vignettes/%.html: ../simcross/vignettes/%.Rmd assets/vignettes/badges.html assets/vignettes/paste_badges.R
	cd $(<D); \
	R -e "rmarkdown::render('$(<F)')"; \
	mv $(@F) ../../Web/$(@D)/
	cd $(@D);paste_badges.R $(@F)
