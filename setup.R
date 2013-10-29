## Copyright 2013 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## It is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## See <http://www.gnu.org/licenses/>

###############################################################################
## INITIAL SETUP
###############################################################################

message("install packages")

## fetch Bioconductor installation script
source("http://bioconductor.org/biocLite.R")
## load package names
source("include/packages.R")
source("include/path.R")

## test R version
if (R.version$major < 3) {
  stop("R >= 3.0.0 is needed!")
}

## cleaver needs Bioconductor >= 2.13
if (BiocInstaller::biocVersion() < "2.13") {
  stop("Bioconductor >= 2.13 needed!")
}

## install needed packages from Bioconductor/CRAN
packages <- c("BRAIN", "cleaver", packages)
biocLite(packages)

###############################################################################
## END OF INITIAL SETUP
###############################################################################


###############################################################################
## FETCH EXPERIMENTAL DATA
###############################################################################

message("download data")

## download data from figshare
figshareUrl <- paste0("http://files.figshare.com/",
  c("1207617/referencetable.csv",
    "1207618/pcr.fasta",
    "1207619/drosophila_melanogaster_reference_proteome.fasta",
    "1207620/figure7.zip",
    "1207625/smm.zip"))
destFile <- file.path(dataPath, basename(figshareUrl))

for (i in seq(along=figshareUrl)) {
  download.file(url=figshareUrl[i], destfile=destFile[i])

  if (grepl(pattern="\\.zip$", x=destFile[i])) {
    ## unzip zip files
    unzip(zipfile=destFile[i], exdir=dataPath, overwrite=TRUE)
    ## remove zip files
    unlink(destFile[i])
  }
}

###############################################################################
## END OF FETCH EXPERIMENTAL DATA
###############################################################################

