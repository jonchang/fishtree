all: readme

readme:
	Rscript -e 'rmarkdown::render("README.Rmd")'
