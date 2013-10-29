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
## SETUP
###############################################################################

## show warnings when they appear
options(warn=1)

## load packages
library("MALDIquant")
library("MALDIquantForeign")

## load include files
source(file.path("include", "path.R"))
source(file.path("include", "plot-functions.R"))

###############################################################################
## END OF SETUP
###############################################################################


###############################################################################
## DATA IMPORT
###############################################################################

## import spectra
message("data import")
## (suppress HPC warnings. See readBrukerFlexData:::.hpc for details.)
rawSpectra <- suppressWarnings(import(file.path(dataPath, "figure7"),
                                      verbose=TRUE))
names(rawSpectra) <- sapply(rawSpectra, function(x)metaData(x)$name)

###############################################################################
## END OF DATA IMPORT
###############################################################################


###############################################################################
## PREPROCESSING
###############################################################################

## variance stabilization
message("variance stabilization")
transformedSpectra <- transformIntensity(rawSpectra, method="sqrt")

## smoothing with a 7-point moving average filter
message("smoothing")
smoothedSpectra <- smoothIntensity(transformedSpectra, method="MovingAverage",
                                   halfWindowSize=3)

## baseline correction
message("baseline correction")
baselineCorrectedSpectra <- list(
  "IPP.L5" = removeBaseline(smoothedSpectra[["IPP.L5"]], method="TopHat",
                            halfWindowSize=1000),
  "SMM.K4" = removeBaseline(smoothedSpectra[["SMM.K4"]], method="TopHat",
                            halfWindowSize=20))

## intensity normalization/calibration using the TIC
message("intensity calibration")
normalizedSpectra <- calibrateIntensity(baselineCorrectedSpectra, method="TIC")

###############################################################################
## END OF PREPROCESSING
###############################################################################


###############################################################################
## FIGURES
###############################################################################

message("plotting")
## 1-column figure could be 85mm: 85mm/25.4
## 2-column figure could be 170mm: 170mm/25.4
oneColumnInch <- 75/25.4
twoColumnInch <- 150/25.4

## 3.35 inches is very small; multiply by a scale factor
oneColumnWidth <- oneColumnInch * 2
twoColumnWidth <- twoColumnInch * 2

## set default EPS options
setEPS()

## set colors
col <- c(IPP.L5="#008040", SMM.K4="#008040")
for (i in seq(along=normalizedSpectra)) {
  metaData(normalizedSpectra[[i]]) <-
    modifyList(metaData(normalizedSpectra[[i]]), list(col=col[i]))
}

## plot IPP
postscript(file.path(plotPath, "IPP.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plotSpectra(normalizedSpectra["IPP.L5"], xlim=c(1450, 16000))
dev.off()

## plot IPP vs SMM 1600-4020 m/z
postscript(file.path(plotPath, "IPPvsSMM_1600-4020.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plotSpectra(normalizedSpectra, xlim=c(1600, 4020))
dev.off()

###############################################################################
## END OF FIGURES
###############################################################################

