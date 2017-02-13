
render:
	Rscript -e "rmarkdown::render_site()"

publish:
	mkdir -p /virdir/Scratch/wwwexport/BBMRIomics/
	cp -R _site/* /virdir/Scratch/wwwexport/BBMRIomics/

clean: 
	Rscript -e "rmarkdown::clean_site()"
