all: doc vignettes data
.PHONY: doc vignettes data

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

# build package documentation
doc:
	R -e 'devtools::document()'

vignettes: inst/doc/simcross.html

inst/doc/simcross.html: vignettes/simcross.Rmd
	cd $(@D);R ${R_OPTS} -e 'library(knitr);knit2html("../../$<")'

data: data/AILped.RData

data/AILped.RData: inst/scripts/create_AILped_dataset.R
	cd $(<D);R CMD BATCH ${R_OPTS} $(<F)
