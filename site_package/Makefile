render:
	Rscript -e "rmarkdown::render_site()"

publish:
	mkdir -p /virdir/Scratch/wwwexport/BBMRIomics/
	cp -R _site/* /virdir/Scratch/wwwexport/BBMRIomics/

clean: 
	Rscript -e "rmarkdown::clean_site()"

newpage:
	cp _template.Rmd $(file).Rmd

apps:
	mkdir -p /srv/shiny-server/BBMRIomics/DNAm/
	ln -snf /opt/shiny-server/BBMRIomics/DNAm/* /srv/shiny-server/BBMRIomics/DNAm/
