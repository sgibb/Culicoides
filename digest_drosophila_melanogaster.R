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

## load packages
library("BRAIN")
library("cleaver")

## load include files
source(file.path("include", "path.R"))
source(file.path("include", "hdpr-functions.R"))
source(file.path("include", "plot-functions.R"))

###############################################################################
## END OF SETUP
###############################################################################


###############################################################################
## DATA IMPORT
###############################################################################

## load Drosophila Melanogaster Reference Proteome
message("protein data import")
proteins <- Biostrings::readAAStringSet(
  file.path(dataPath, "drosophila_melanogaster-reference_proteome.fasta"),
  format="fasta")

###############################################################################
## END OF DATA IMPORT
###############################################################################


###############################################################################
## PREPROCESSING
###############################################################################

## in silico cleavage
message("in silico cleavage of reference proteins (be patient)")
peptides <- cleaver::cleave(proteins, enzym="trypsin", missedCleavages=0:1)
peptides <- as.character(unlist(peptides), use.names=FALSE)
# length(peptides) == 2657056

## some peptides contain selenocysteine (U; n=20) which is not supported by
## BRAIN and some peptides contain an unknown amino acid (X; 478)
peptides <- peptides[!grepl(pattern="U|X", x=peptides)]
# length(peptides) == 2656558

## calculate mass values for peptide sequences
drosophilaMelanogasterMass <- sapply(peptides, function(x) {
  return(BRAIN::calculateMonoisotopicMass(BRAIN::getAtomsFromSeq(x)))
})

## use only mass in 700:4000 Da
drosophilaMelanogasterMass <- drosophilaMelanogasterMass[
  drosophilaMelanogasterMass >= 700 &
  drosophilaMelanogasterMass <= 4000]
# length(drosophilaMelanogasterMass) == 1595674

## find HDPR regression line
hdprDrosophilaMelanogasterCoef <- hdprCoef(drosophilaMelanogasterMass)
message("HDPR-line: y=", round(hdprDrosophilaMelanogasterCoef[2], 5), "*x+",
        round(hdprDrosophilaMelanogasterCoef[1], 5))

diffReg <- abs((hdprDrosophilaMelanogasterCoef[1]+
            hdprDrosophilaMelanogasterCoef[2]*drosophilaMelanogasterMass) -
            decimalPlaceHdpr(drosophilaMelanogasterMass))

print(quantile(diffReg, probs=c(0.5, 0.9, 0.95, 0.97, 0.98, 0.99)))

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

## plot half decimal place rule
postscript(file.path(plotPath, "hdpr_drosophila_melanogaster.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plotHdpr(drosophilaMelanogasterMass, hdprDrosophilaMelanogasterCoef, tolerance=0.2)
dev.off()

## plot cdf
postscript(file.path(plotPath, "cdf_drosophila_melanogaster.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
plot(ecdf(diffReg))
abline(h=0.98, v=0.2, col=2, lty=3)
dev.off()

###############################################################################
## END OF FIGURES
###############################################################################

