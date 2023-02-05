# BBMRIomics

BBMRIomics is an R package that facilitates BBMRI-omics downstream
analysis that is availabe at the [BBMRIomics workspace](TBA) and previously the so-called
[BBMRI BIOS virtual machine](http://www.bbmriwiki.nl/wiki/BIOS_VirtualMachine)
running at SURF Research Cloud.

For an introduction and examples, visit
[BBMRIomics](http://bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/).

See also the [TODO-list](TODO.md).

## Features

* Easy access to BBMRI-omics data; RNAseq, DNAm, metabolomics and
  sequenced- and imputed-genotypes.
* Contains a collection of preprocessed data sets in so-called
  [`SummarizedExperiments`](http://bioconductor.org/packages/SummarizedExperiment/).
* These `SummarizedExperiments` contain both the actual data,
  i.e. counts, beta-values and M-values, as well as metadata on feature- 
  and sample-level.
* Multiomics integrated analysis is simplified by introducing
  an anonymized cross omic identifier.
* Furthermore, links to BBMRI metabolomics or GoNL are provided as
  well.
* Interaction with BBMRIomics underlying metadatabase is provided
  through a number of functions, such as the `getSQLview`-function.

# For users
Here is a step-by-step manual to start working on BBMRIomics datasets which is provisioned at [SURF's Research Cloud](https://www.surf.nl/en/surf-research-cloud-collaboration-portal-for-research). A Code of Conduct (CoC) document needs to be signed first and you can contact [Rick Jansen](https://psychiatryamsterdam.nl/researcher/jansen-rick) for this. If you have any question, please contact Leon Mei at LUMC.

A video instruction to show how to access and use SURF's Research Cloud from [a researcher's perspective](https://www.youtube.com/watch?v=_rdK3W1AWvA)

## Use of SURF research access management (SRAM)
You can log into all web portals described below with your own institutional account thanks to the support of SRAM except for the access to the actual workspace where you will be using a Time based password. 

## Step 1: Request personal budget
* First, please log into https://servicedesk.surf.nl/ and fill in the form of "Small NWO request (EINF)". 
* If you are a PhD student, the "Signing Authority" should be your supervisor's contact details.
* "Resources" should be "Research Cloud - HPC Cloud"
* Please mention BBMRIomics resource will be used in the "Description" field
* You can also request up to 160 consultancy hours for "Research Cloud - HPC Cloud: how many consultancy hours (default 8 hours)" if you think you will need additional SURF support on running the Research Cloud services 
* CPU and GPU core hours can be filled in based on your specific needs. Normally we only need CPU and you could put the maximum number of hours of "50000" in this application. 

## Step 2: Configuration of SRAM BBMRIomics CO
* After you have signed CoC, you will get an invitation to this CO, log into https://sram.surf.nl/ using your institutional credential (using "Login via SURFconnext") to verify you are a member of SRC-BBMRI CO.
* Alternatively, user can request access to "BBMRIomics CO" through this [link](https://sram.surf.nl/registration?collaboration=8051068c-59be-4580-b67b-210cfd02dd89)
* Optionally if you will need to access the workspace's Linux terminal, you can set up your [SSH key in your user profile section](https://servicedesk.surf.nl/wiki/display/WIKI/Log+in+to+your+workspace#Logintoyourworkspace-AccessaworkspacewithSSH).

## Step 3: Configuration of Research Drive (RD)
* You will get an invitation email to join Research Drive. After that, log into https://researchdrive.surfsara.nl using your institutional credential. Choose "Login via SURFconnext". You may get an error message saying "Oops! Your SURF Research Drive account with email address `<xxxx>` and account `<yyyyyy>-lumcnet@lumc.nl` could not be found. It seems that you do not have permission to use SURF Research Drive. Please contact the IT department of your institution or data steward of your project to request access to SURF Research Drive.". In this case, please forward this message to Leon Mei to so that he can properly share the datasets.
* You should see a shared folder called "RSC_BIOS". There are two sub-directories "RP3_data" and "RP3_analysis" that contains the prepared BBMRIomics data for downstream analysis. They are read-only and you need to store the analysis results in your own personal space on RD. 
* **Please be aware that there is no backup support of your data stored at Research Drive, so make sure you arrange proper backup of your scripts and data yourself.**

## Step 4: Installation and configuration of authenticator on your mobile phone
* Research Cloud is using a so called [Time-based One-time Password authentication (TOTP) for authentication](https://servicedesk.surf.nl/wiki/display/WIKI/Log+in+to+your+workspace#Logintoyourworkspace-WorkspaceAccesswithTOTP). 
* You need to install an authenticator app on your mobile phone, e.g., from [Google appstore](https://play.google.com/store/apps/details?id=com.google.android.apps.authenticator2) or [Apple appstore](https://apps.apple.com/us/app/google-authenticator/id388497605).
* In this authenticator app, you will see a username and a numeric password (that is changing every minute). This username and password are used to log into the R-studio environment. 

## Step 5: Configuration of your account at Research Cloud
* Log into https://portal.live.surfresearchcloud.nl/ using your institutional credential
* If you don't have budget, you can click "Request" button in "Request a new wallet" to apply for more resources using your own Personal EINF Budget.
* Go to "Profile" tab using the navigation bar at the top of the page.
  * Set up time based password in the "Profile" section by clicking the "more option" button. [more details](https://servicedesk.surf.nl/wiki/display/WIKI/Log+in+to+your+workspace#Logintoyourworkspace-WorkspaceAccesswithTOTP).
  * Link Researchdrive in the "Collaborative organisations" section by by clicking "more option" button. [more details](https://servicedesk.surf.nl/wiki/display/WIKI/Connect+Research+Drive)  
* Once your request to budget/new wallet gets approved, you can request your user specific persistent storage and fixed IP address. This persistent storage and IP resource can be used when you create a new workspace.
* Additional documentations about Research Cloud: https://servicedesk.surf.nl/wiki/display/WIKI/Research+Cloud+Documentation

## Step 6: Launch or access a BBMRIomics workspace
* In the main dashboard of Research Cloud, click "Create new workspace" and follow the instructions to configure and start a new workspace.
  * **About the expiration date, please try to keep a reasonable timeline. If the workspace is running without doing any jobs, it is a waste of the Cloud resource**
  * **When you are not using a workspace, you can also pause the workspace so that it won't consume resources.** 

More information can be found at https://servicedesk.surf.nl/wiki/display/WIKI/Start+a+simple+workspace and https://servicedesk.surf.nl/wiki/display/WIKI/Workspaces

* As a member of the same CO, you can also access a running workspace created by another member in the same CO provided the owner agrees.
* After you log into the Rstudio UI, you can run the following commands to verify if you have proper access to the BBMRIomics data.
```{r}
library(BBMRIomics)
data(package="BBMRIomics")
bbmri.data("rp3_rp4_meta")
head(person)
```

## Step 7: Manage personal (R) packages
* Some basic R packages are available in the workspace. For other common packages, you could request their installation via Leon Mei or Davy Cats.
* You can also install your own libraries in the running workspace in the persistent storage.
* All data/scripts stored in your researchdrive folder are persistent (although no backup support).
* More user documents about the BBMRIomics package can be found at https://bbmri-nl.github.io/BBMRIomics/usage.html

## (optional) Step 8: ssh access
For more advanced users who would like to access the terminal of a workspace. You can `ssh <user_name>@<IP>` to the workspace where <user_name> is the same as defined in your authenticator app and the <IP> address is shown in the "Workspace Details". You can "sudo su" to become a sudoer on the workspace in order to install necessary programs. 

# For developers 

The **BBMRIomics** git repo contains three directories: 

1. `BBMRIomics`: containing the actual **BBMRIomics** package
2. `site_package`: containing Rmarkdown code for generation of
   [BBMRIomics](bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/index.html)
   vignettes collected as a web-site.
3. `couchapp`: a directory containing the couchdb metadatabase
   JSON-schema and some scripts to use
   [couchapp](https://github.com/couchapp/couchapp) for interaction
   with the metadatabase (DEPRECATED)

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

After installation the preprocessed datasets need to be linked for easy loading
the data i.e., using `data(<dataset>)`. The `inst/configure` or `configure` 
directory contains a `Makefile` to accomplish this.

```bash
sudo -i make -f /opt/R/microsoft-r/3.3/lib64/R/library/BBMRIomics/configure/Makefile
```

## Guidelines for contributing to the web-site ##

We welcome contributions to the
[BBMRIomics](bios-vm.bbmrirp3-lumc.surf-hosted.nl/BBMRIomics/index.html)
web-site, for example, with use-cases. The BBMRIomics web-site is
completely build using R. Specifically,
[Rmarkdown](http://rmarkdown.rstudio.com/) is used to generate the
complete web-site from a few Rmarkdown -files. This Rmarkdown
[section](http://rmarkdown.rstudio.com/rmarkdown_websites.html) shows
how easy it is to generate a web-site using rmarkdown.

### Adding a use-case ###

Modify `_site.yml` by adding a new use-case under the use-case items, e.g. 

- text: "eqtl analysis"
  href: eqtl.html
  
*Beware: spacing does matter!*

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

and if it runs successfully

```bash
sudo make publish
```

### Creating datasets ###

The `developers`-section of the web-site presents reproducible
vignettes to construct the datasets available via the
*BBMRIomics*-package. These vignettes use a different setup-mechanism,
compared to the example and use-case vignettes, to ensure they are
only build once or specifically if the code-chunk option `eval` is set
to `TRUE` by a developer.



