all: doc data
.PHONY: doc data

# R_OPTS: --vanilla without --no-environ
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

# build package documentation
doc:
	R -e 'devtools::document()'

data: data/AILped.RData

data/AILped.RData: inst/scripts/create_AILped_dataset.R
	cd $(<D);R CMD BATCH ${R_OPTS} $(<F)
