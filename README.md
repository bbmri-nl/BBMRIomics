# BBMRIomics

**BBMRIomics** is an R package that facilitates BBMRI-omics downstream
analysis using R on surfSARA HCP virtual machines.

For an introduction and examples, visit the [BBMRIomics](http://bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/).

## Features

* Easy access to BBMRI-omics data; RNAseq, DNAm, metabolomics and sequenced- and imputed-genotypes.
* Contains a collection of preprocessed data sets so-called [`SummarizedExperiments`](http://bioconductor.org/packages/SummarizedExperiment/).
* These `SummarizedExperiments` contain both the actual data, i.e. counts, beta-values, as well as metadata on feature- and sample-level.
* Across omic-level integrated analysis is simplified by introducing an anonymized cross omic identifier.
* Furthermore, links to BBMRI metabolomics or GoNL are provided as well.
* Access to all BBMRI BRAINSHAKE metabolomics data is provided by an R-interface to `molgenis` RESTfull-API.
* Interaction with BBMRIomics underlying metadatabase is provided through a `view`-function.

## User requirements

Currently, access to the metadatabase and `molgenis` RESTfull-API
requires additional accounts (contact Leon Mei). To avoid repeately
typing of usernames and passwords **BBMRIomics** uses a configuration
file (`~/.bbmriomics`) which should be place in your home-directory on
the VM and containing the account of the metadatabase (`usrpwdrp3`) or
molgenis metabolomics (`usrpwdrp4`) and optionally the location of
your grid_proxy:

usrpwdrp3: '<usrname:password>'
usrpwdrp4: '<usrname:password>'
proxy: /tmp/grid_proxy

One loading the **BBMRIomics** in R

```r
library(BBMRIomics)
```

the configuration file will be loaded and your account setting will be
store for easy use in this R session. For example, connecting to the `molgenis` metabolomics database:

```r
molgenis.connect(usrpwd=RP4_DB_USRPWD, url=RP4_DB)
```

## Development notes

This git repo contains the **BBMRIomics** package, `site_package` code
to generate the
[BBMRIomics](bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/index.html)
web-site and a directory containing the `couchdbapp` for easy
interaction with the **BBMRIomics** underlying metadatabase.

### Installation **BBMRIomics**-package ###

Installation of the **BBMRIomics** requires sudo-rights. Furthermore,
use `sudo -i` to load the environmental variable `R_SITE_LIB` which
should be `/opt/...` to install the package out-side the VM's storage
space. Either use R CMD build, R INSTALL or from within R 

```r
library(devtools)
setwd("/BBMRIomics")
install()
```

After installation views need to be updated and the preprocessed
datasets need linked for easy loading (using `data()`). The tools
directory contains a `Makefile` to accomplish this.

### Guidelines for contributing to the web-site ###

We welcome contributions to the
[BBMRIomics](bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/index.html)
web-site, for example, with use-cases. The **BBMRIomics** web-site is
fully build using R. Especially, we use
[Rmarkdown](http://rmarkdown.rstudio.com/) to generate the complete
web-site. The
[section](http://rmarkdown.rstudio.com/rmarkdown_websites.html) show
how easy it is to generate from a collection of **rmarkdown**-files a
website. Basically the only requirement is a `_site.yml`-file. Of
course you should be familiar with **rmarkdown**. 

#### Adding a use-case ####

Modify `_site.yml` by adding a new use-case under the use-case items, e.g. 
- text:  "eqtl analysis"
  href: eqtl.html
  
*beware spacing does matter!*

Now from within the `site_package` directory generate a new template using 

```bash
make newpage file=eqtl
```

This will generate `eqtl.Rmd` for you with the same settings as the
other Rmd-files in this folder. Provide your example with
documentation and use `cache=TRUE` as an argument for the R-code
chunk (maybe check one of the other Rmd-documents). 

*Check if all code chunks run without errors, optionally add a library to `_setup.Rmd`*. Now run

```bash
make render     
```

and if successfully run

```bash
sudo make publish
```
, which requires sudo rights.


### Interaction with the metadatabase ###


## TODO's ##


