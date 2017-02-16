render:
	Rscript -e "rmarkdown::render_site()"

publish:
	#mkdir -p /virdir/Scratch/wwwexport/BBMRIomics/
	#cp -R _site/* /virdir/Scratch/wwwexport/BBMRIomics/
	cp -R _site /virdir/Scratch/mviterson

clean: 
	Rscript -e "rmarkdown::clean_site()"

newpage:
	cp _template.Rmd $(file).Rmd
