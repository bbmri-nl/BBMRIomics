This README describes additional configuration and installation of
RStudio server/shiny server and R packages on the surfSARA HPC Virtual
Machines

last editted by: Maarten van Iterson
date: 15 March 2017

# Install R packages #

clean Ubuntu installation requires a few additional libraries.

required for `openssl` which is required for R package `devtools`:

```bash
> sudo apt-get install libssl-dev
```

required for R package `XML`

```bash
> sudo apt-get install libxml2-dev 
```
  
when installing using devtools first run:

```r
httr::set_config(httr::config( ssl_verifypeer = 0L))
devtools::install_github(<github_package>)
```

# Additional R/RStudio Server configuration #

Overwrite default R library path:

add to `/etc/rstudio/rsession.conf`
```bash
r-libs-user=/opt/R/microsoft-r/3.3/lib64/R/library
```

add to `R.home()/etc/Rprofile.site`

```bash
export R_LIBS_SITE=/opt/R/microsoft-r/3.3/lib64/R/library

.First <- function() {
cat("\n   setting MKL threads to 2!\n\n")
RevoUtilsMath::setMKLthreads(2)
##optionally set user memory limit
ulimit::memory_limit(10240) ##in MiB = 10Gb
}

#do not ask for saving and prevent load .RData/.Rhistory

utils::assignInNamespace("q", 
   function(save = "no", status = 0, runLast = TRUE) {
     .Internal(quit(save, status, runLast))
   }, "base")    

utils::assignInNamespace("quit", 
   function(save = "no", status = 0, runLast = TRUE) {
     .Internal(quit(save, status, runLast))
   }, "base")
```


# Shiny Server configuration #

TODO



