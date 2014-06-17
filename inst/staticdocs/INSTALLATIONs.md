## 1. R requirement

R (http://www.r-project.org) is a language and environment for statistical computing and graphics. We assume R (version 3.1.0) has been installed in your local machine. It can be installed following quick instructions below for different platforms (Windows, Mac, and Linux).

* Quick link for `Windows`: [Download R for Windows](http://www.stats.bris.ac.uk/R/bin/windows/base/R-3.1.0-win.exe).
* Quick link for `Mac`: [Download R for Mac OS X 10.6 (Snow Leopard or higher)](http://cran.r-project.org/bin/macosx/R-3.1.0-snowleopard.pkg).

* Below are `shell command lines in Terminal` (for `Linux`):

Assume you have a `ROOT (sudo)` privilege:
    
    sudo su
    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.1.0.tar.gz
    tar xvfz R-3.1.0.tar.gz
    cd R-3.1.0
    ./configure
    make
    make check
    make install
    R # start R

Assume you do not have a ROOT privilege and want R installation under your home directory (below `/home/hfang` should be replaced with yours):

    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.1.0.tar.gz
    tar xvfz R-3.1.0.tar.gz
    cd R-3.1.0
    ./configure --prefix=/home/hfang/R-3.1.0
    make
    make check
    make install
    /home/hfang/R-3.1.0/bin/R # start R

## 2. Installation of the package

Notes: below are `R command lines (NOT shell command lines in Terminal)`.

First, install dependant/imported/suggested packages:

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("hexbin","ape","supraHex","graph","Rgraphviz","igraph","Biobase","limma","survival","foreach","doMC"))

Second, install the package `dcGOR` under [stable release version](http://cran.r-project.org/package=dcGOR):

    install.packages("dcGOR",repos="http://cran.r-project.org",type="source")

Notes: to install the package `dcGOR` under [latest development version](https://github.com/hfang-bristol/dcGOR), please run:

    source("http://bioconductor.org/biocLite.R")
    biocLite("devtools")
    library(devtools)
    install_github(c("dcGOR"), username="hfang-bristol")

## 3. Workflow of the package

It provides a brief overview of how the package operates and what you expect to get from it.

<a href="javascript:newWin('dcGOR_workflow.png', 'dcGOR_workflow.png', '1200', '1200')" title="Click to enlarge"><img style="max-width:95%;border:1px solid #0000FF;box-shadow:5px 5px 2px #C0C0FF;" src='dcGOR_workflow.png', width="800px" /></a>
