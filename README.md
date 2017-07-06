# BBMRIomics

last editted by: Maarten van Iterson

date: 30 June 2017


BBMRIomics is an R package that facilitates BBMRI-omics downstream
analysis that is availabe at [BBMRI BIOS virtual machine](http://www.bbmriwiki.nl/wiki/BIOS_VirtualMachine)
running at surfSARA HPC Cloud.

For an introduction and examples, visit
[BBMRIomics](http://bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/).

## Features

* Easy access to BBMRI-omics data; RNAseq, DNAm, metabolomics and
  sequenced- and imputed-genotypes.
* Contains a collection of preprocessed data sets so-called
  [`SummarizedExperiments`](http://bioconductor.org/packages/SummarizedExperiment/).
* These `SummarizedExperiments` contain both the actual data,
  i.e. counts, beta-values, as well as metadata on feature- and
  sample-level.
* Multiomics integrated analysis is simplified by introducing
  an anonymized cross omic identifier.
* Furthermore, links to BBMRI metabolomics or GoNL are provided as
  well.
* Access to all BBMRI BRAINSHAKE metabolomics data is provided by an
  R-interface to molgenis RESTful-API.
* Interaction with BBMRIomics underlying metadatabase is provided
  through a `view`-function.

## User requirements

Currently, access to the metadatabase and molgenis RESTful-API
requires additional accounts (please check [BBMRI-BIOS wiki](http://www.bbmriwiki.nl/wiki/BIOS_VirtualMachine#BIOSVMAccess)
for instructions). To avoid repeately
typing of usernames and passwords BBMRIomics uses a configuration file
(`~/.bbmriomics`) which should be placed in your home-directory on the
VM and containing the account of the metadatabase (`usrpwdrp3`) or
molgenis metabolomics database (`usrpwdrp4`) and optionally the
location of your grid_proxy if you want to use the `SRM2VM`
file-download function from within R. Finally, the configuration file
(`~/.bbmriomics`) should look like this:

```{bash}
usrpwdrp3: '<usrname:password>'
usrpwdrp4: '<usrname:password>'
proxy: /tmp/<grid_proxy>
```

When loading the **BBMRIomics**- package in R

```r
library(BBMRIomics)
```

the configuration file will be loaded and your account settings will be
stored for easy use in the current R session. For example, connecting to the
molgenis metabolomics database using the stored `RP4_DB_USRPWD` read from 
the configuration file, is done like this:

```r
molgenis.connect(usrpwd=RP4_DB_USRPWD, url=RP4_DB)
```

> NOTE: urls and data-directories are stored in the package
> subdirectory `inst/configure`. If urls or data directories change in
> the future this file should be modified.

# For developers 

The **BBMRIomics** git repo contains three directories: 

1. `BBMRIomics`: containing the actual **BBMRIomics** package
2. `site_package`: containing Rmarkdown code for generation of
   [BBMRIomics](bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/index.html)
   vignettes collected as a web-site.
3. `couchapp`: a directory containing the couchdb metadatabase
   JSON-schema and some scripts to use
   [couchapp](https://github.com/couchapp/couchapp) for interaction
   with the metadatabase (IN DEVELOPMENT)

The following sections describe the use and maintenance of each in
detail.

## Installation of the **BBMRIomics**-package ##

Installation of the **BBMRIomics**-package on BIOS VMs requires
sudo-rights. Furthermore, use `sudo -i` to load the environmental
variable `R_SITE_LIB`, which is for the current R installation
`/opt/R/microsoft-r/3.3/lib64/R/library`, to install the package
persistently, i.e., out-side possibly temporarily VM storage-space.

A clean installation of Ubuntu requires the installation of some
additional headers and libraries. See
[prerequisites](BBMRIomics/inst/configure/prerequisites.md) for additional
instructions.

In general R libraries are installed using: 

```{bash}
sudo -i R
```
and then from within R:

```{r}
library(BiocInstaller)
biocLite(<pkgName>)
```

This allows installation of packages from CRAN, BioConductor or
github.

### Installation from git clone locally ###

Either use from the command-line: 

```{bash}
R CMD build BBMRIomics
sudo -i R CMD INSTALL BBMRIomics_x.y.z.gz
```

or more easily from within R, started with `sudo -i R`:

```{r}
library(devtools)
setwd("/BBMRIomics")
install()
```

### Installation from github directly ###

```{r}
library(devtools)
install_github("bbmri-nl/BBMRIomics", subdir="BBMRIomics")
```

After installation the views need to be updated and the preprocessed
datasets need to be linked for easy loading the data i.e., using
`data(<dataset>)`. The `inst/configure` or `configure` directory
contains a `Makefile` to accomplish this.

```bash
sudo -i make -f /opt/R/microsoft-r/3.3/lib64/R/library/BBMRIomics/configure/Makefile USRPWD='<usrpwdrp3>'
```

## Guidelines for contributing to the web-site ##

We welcome contributions to the
[BBMRIomics](bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/index.html)
web-site, for example, with use-cases. The BBMRIomics web-site is
completely build using R. Especially,
[Rmarkdown](http://rmarkdown.rstudio.com/) is used to generate the
complete web-site from a few Rmarkdown -files. This Rmarkdown
[section](http://rmarkdown.rstudio.com/rmarkdown_websites.html) shows
how easy it is to generate a web-site using rmarkdown.

### Adding a use-case ###

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

Check if all code chunks run without errors, optionally add a library
to `_setup.Rmd`. Now run

```bash
make render     
```

and if successfully run

```bash
sudo make publish
```

## Interaction with the metadatabase ##







