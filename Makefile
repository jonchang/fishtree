all: readme

readme:
	Rscript -e 'devtools::build_readme()'
