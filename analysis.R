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
source(file.path("include", "packages.R"))
lapply(packages, library, character.only=TRUE)

## load include files
lapply(list.files("include", pattern="^.*\\.R$", full.names=TRUE), source)

###############################################################################
## END OF SETUP
###############################################################################


###############################################################################
## DATA IMPORT
###############################################################################

## load reference table
message("data import")
mapFile <- list.files(path=dataPath, recursive=TRUE,
                      pattern="referencetable.csv")[1]
map <- read.csv(file.path(dataPath, mapFile), sep=";", dec=".",
                stringsAsFactors=FALSE)

## import spectra
## (suppress HPC warnings. See readBrukerFlexData:::.hpc for details.)
rawSpectra <- suppressWarnings(import(file.path(dataPath, unique(map$folder)),
                                      verbose=TRUE))

## rewrite some metaData
rawSpectra <- lapply(rawSpectra, function(x) {
  ## identify spectrum by spot & folder
  spot <- metaData(x)$spot
  file <- metaData(x)$file

  lSpot <- map$spot == spot
  lFile <- unlist(lapply(map$folder, grepl, x=file, fixed=TRUE))

  i <- which(lSpot & lFile)

  metaData(x)$name <- map[i, "pcr"]
  metaData(x)$label <- map[i, "label"]
  metaData(x)$species <- substring(map[i, "label"], 1, 5)
  metaData(x)$group <- switch(metaData(x)$species,
                              "C_Pun" =,
                              "C_New" =,
                              "C_Pul" = { "G_Pul" },
                              "C_Nub" = { "G_Nub" },
                              "C_Dew" =,
                              "C_Sco" =,
                              "C_Obs" = { "G_Obs" },
                              { "G_UNK" })
  metaData(x)$col <- switch(metaData(x)$species,
                              "C_Pun" = "#0080FF",
                              "C_New" = "#00009B",
                              "C_Pul" = "#51519B",
                              "C_Nub" = "#008040",
                              "C_Dew" = "#FF0000",
                              "C_Sco" = "#B30000",
                              "C_Obs" = "#FFB0AE",
                              { "#000000" })

  metaData(x)$spot <- paste(map[i, "folder"], map[i, "spot"], sep="-")
  metaData(x)$isCalibration <- map[i, "name"] == 0

  return(x)
})

## remove calibration spectra
isCalibration <- unlist(lapply(rawSpectra, function(x) {
  metaData(x)$isCalibration
}))
rawSpectra <- rawSpectra[!isCalibration]

## load pcr FASTA data
pcr <- ape::read.dna(file.path(dataPath, "pcr.fasta"), format="fasta")

###############################################################################
## END OF DATA IMPORT
###############################################################################


###############################################################################
## PREPROCESSING
###############################################################################

## check raw data
message("check raw data ")
isAnyEmpty <- any(sapply(rawSpectra, isEmpty))
isAnyIrregular <- any(!sapply(rawSpectra, isRegular))
## spectra measured at 2012-11-20 contain 114200 datapoints
## spectra measured at 2013-04-05 contain 114222 datapoints
## this difference doesn't matter
isAnyLengthDifferent <- !all(
  (length(rawSpectra[[1]]) == sapply(rawSpectra, length)) |
  (length(rawSpectra[[length(rawSpectra)]]) == sapply(rawSpectra, length)))

areAllSpectraValid <- !(isAnyEmpty | isAnyIrregular | isAnyLengthDifferent)
stopifnot(areAllSpectraValid)

## variance stabilization
message("variance stabilization")
transformedSpectra <- transformIntensity(rawSpectra, method="sqrt")

## smoothing with a 7-point moving average filter
message("smoothing")
smoothedSpectra <- smoothIntensity(transformedSpectra, method="MovingAverage",
                                   halfWindowSize=3)

## baseline correction
message("baseline correction")
baselineCorrectedSpectra <- removeBaseline(smoothedSpectra,
                                           method="TopHat",
                                           halfWindowSize=20)

## intensity normalization/calibration using the TIC
message("intensity calibration")
normalizedSpectra <- calibrateIntensity(baselineCorrectedSpectra, method="TIC")

## peak detection; noise estimation by Friedman's SuperSmoother
message("peak detection")
peaks <- detectPeaks(normalizedSpectra, halfWindowSize=20, SNR=2,
                     method="SuperSmoother")

## alignment/recalibration of spectra using a quadratic warping function
message("peak warping")
warpingFunctions <- determineWarpingFunctions(peaks, minFrequency=0.85,
                                              tolerance=200e-6,
                                              method="quadratic")

## apply warping functions
warpedSpectra <- warpMassSpectra(normalizedSpectra, warpingFunctions)
warpedPeaks <- warpMassPeaks(peaks, warpingFunctions)

## run monoisotopic peak detection
message("monoisotopic peak detection")
# load averagine reference table
load(file.path("include", "averagineTable.RData"))

monoisotopicPeaks <- monoisotopic(warpedPeaks, chargeState=1, tolerance=50e-6,
                                  intensityTolerance=0.3, SNR=2,
                                  referenceTable=averagineTable)

## unique mass (needed for HDPR plotting)
uniqueMass <- unique(unlist(lapply(monoisotopicPeaks, mass)))

## remove proton (H+)
uniqueMass <- uniqueMass - 1.0078250321

## find HDPR regression line
hdprCo <- hdprCoef(uniqueMass)
message("HDPR-line: y=", round(hdprCo[2], 5), "*x+", round(hdprCo[1], 5))

## filter using HDPR
hdprFilteredPeaks <- lapply(monoisotopicPeaks, function(p) {
  return(p[hdprFilter(mass(p), co=hdprCo, tolerance=0.2)])
})

## peak binning
message("peak binning")
binnedPeaks <- binPeaks(hdprFilteredPeaks, tolerance=200e-6, method="strict")

## create labels
labels <- unlist(lapply(binnedPeaks, function(x) {
  return(metaData(x)$label)
}))

## filter peaks, only peaks occurring in 2 of 3 technical replicates will
## survive
message("peak filtering")
filteredPeaks <- filterPeaks(binnedPeaks, minFrequency=2/3, labels=labels)

## merge technical replicates
message("merge technical replicates")
mergedPeaks <- mergeMassPeaks(filteredPeaks, labels=labels)

## update labels (length of mergedPeaks is 1/3 of binnedPeaks)
labels <- unlist(lapply(mergedPeaks, function(x) {
  return(metaData(x)$label)
}))
speciesLabels <- unlist(lapply(mergedPeaks, function(x) {
  return(metaData(x)$species)
}))
groupLabels <- unlist(lapply(mergedPeaks, function(x) {
  return(metaData(x)$group)
}))

## create peak matrix
message("peak matrix")
peakMatrix <- intensityMatrix(mergedPeaks)
rownames(peakMatrix) <- labels

## create binary peak matrix
binaryPeakMatrix <- ifelse(is.na(peakMatrix), 0, 1)

## save (binary) peak matrix
write.table(peakMatrix, file=file.path(csvPath, "peakmatrix.csv"),
            sep=";", dec=".", col.names=NA)
write.table(binaryPeakMatrix, file=file.path(csvPath, "binarypeakmatrix.csv"),
            sep=";", dec=".", col.names=NA)

## cluster analysis
message("similarity/distance measurement")
similarityMatrix <- proxy::simil(binaryPeakMatrix, method="Dice")
distanceMatrix <- proxy::as.dist(similarityMatrix)

## hierarchical clustering
message("ms clustering & bootstrap (be patient)")
msHclust <- createMsTree(binaryPeakMatrix)
msBoot <- bootstrap(binaryPeakMatrix, createMsTree, n=1e3)

write.table(as.matrix(similarityMatrix),
            file=file.path(csvPath, "similaritymatrix.csv"), sep=";", dec=".",
            col.names=NA)
write.table(as.matrix(distanceMatrix),
            file=file.path(csvPath, "distancematrix.csv"), sep=";", dec=".",
            col.names=NA)
# pcr data
message("pcr clustering & bootstrap (be patient)")
pcr <- pcr[-1, ] ## exclude Forci (outgroup)
pcrHclust <- createPcrTree(pcr)
pcrBoot <- bootstrap(pcr, createPcrTree, n=1e3)

## sda.ranking
message("sda ranking")
rankingMatrix <- peakMatrix
rankingMatrix[is.na(rankingMatrix)] <- 0
## round mass data
colnames(rankingMatrix) <- sprintf("%.3f", as.double(colnames(rankingMatrix)))

speciesRanking <- sda.ranking(rankingMatrix, L=speciesLabels)
groupRanking <- sda.ranking(rankingMatrix, L=groupLabels)

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

## fetch spots
spots <- unlist(lapply(warpedSpectra, function(x) {
  return(metaData(x)$spot)
}))

## plot 121116-I4
i4 <- which(spots == "121116-I4")
i4Peaks <- list(actin=c(1198.71, 1629.82, 1776.79, 1790.89, 1954.07, 2201.06,
                        3151.64),
                "ATP-Synthase alpha subunit (mitochondrial)"=c(1344.73, 1456.81,
                        1553.75, 1610.88),
                "ATP-Synthase beta subunit (mitochondrial)"=c(1435.75, 1601.81, 1617.81,
                        1988.03, 2027.13, 2266.09),
                myosin=c(1159.63, 1206.65, 1586.78, 1819.92, 1941.99,
                        2061.03, 2146.09, 2180.11, 2216.12, 2358.14,
                        2661.31),
                tropomyosin=2218.08)

## WARNING: .which.closest is an internal function of MALDIquant and could
## changed/removed without any information to the user!
i4PeaksIdx <- MALDIquant:::.which.closest(unlist(i4Peaks),
                                          mass(warpedPeaks[[i4]]))
i4Pch <- unlist(mapply(function(p, i)rep(p, length(i)),
                       i=i4Peaks, p=c(4, 2, 6, 1, 0)))

postscript(file.path(plotPath, "I4_peaks-in.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plotSpectra(c(warpedSpectra[[i4]]), main=c(metaData(warpedSpectra[[i4]])$name),
            xlim=c(700, 4000), ylim=c(0, 1.10), add.expr={
            points(warpedPeaks[[i4]][i4PeaksIdx],
                   col="red", pch=i4Pch);
            legend(x=3200, y=0.3, legend=names(i4Peaks), pch=unique(i4Pch),
                   col="red", cex=3/4, bty="n");
            par(new=TRUE, pty="m", plt=c(0.5, 0.95, 0.4, 0.95), cex=0.5, err=-1,
                mgp=c(0, 0.2, 0.2));
            plot(warpedSpectra[[i4]], col=metaData(warpedSpectra[[i4]])$col,
                 main="", sub="", xlim=c(2100, 2280),
                 ylim=c(0, 0.25*max(warpedSpectra[[i4]]@intensity)),
                 xlab="", ylab="", yaxt="n", xaxt="n", bty="n");
            points(warpedPeaks[[i4]][i4PeaksIdx],
                   col="red", pch=i4Pch, cex=1.5);
            axis(1, at=seq(2100, 2280, by=20));
            })
dev.off()

## select representative example spectra
unmergedLabels <- unlist(lapply(filteredPeaks, function(x) {
  return(setNames(metaData(x)$species, metaData(x)$label))
}))

repSpots <- representativeSpots(filteredPeaks, unmergedLabels)

## change order
plotOrder <- c("C_Pul", "C_New", "C_Pun", "C_Nub", "C_Obs", "C_Sco",
               "C_Dew")

repSpots <- repSpots[plotOrder]

## plot representative for each group
postscript(file.path(plotPath, "species.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
plotSpectra(warpedSpectra[repSpots], xlim=c(700, 4000), ylim=c(0, 1))
dev.off()

## plot representative for each group 1820-2900 m/z
topSdaPeaksInRange <- topSdaPeaks(speciesRanking, groupRanking,
                                  xlim=c(1820, 2900), top=c(species=3, group=3))
postscript(file.path(plotPath, "species_1820-2900.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
plotSpectra(warpedSpectra[repSpots],
            speciesMass=topSdaPeaksInRange$species[, plotOrder],
            groupMass=topSdaPeaksInRange$group, xlim=c(1829, 2900),
            ylim=c(0, 0.6))
dev.off()

## plot zoomed pictures for selected spectra
postscript(file.path(plotPath, "zoom.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
## choose the highest entry of 3 species and of all (3) groups
pos <- c(topSdaPeaksInRange$species[1, c("C_Nub", "C_Pul", "C_Sco") ],
         topSdaPeaksInRange$group[1, ])
plotZoomedSpectra(warpedSpectra[repSpots], pos=sort(pos), delta=5,
            speciesMass=topSdaPeaksInRange$species[, plotOrder],
            groupMass=topSdaPeaksInRange$group, ylim=c(0, 0.42))
dev.off()

## colour dendrogram
culicoidesColors <- unlist(lapply(mergedPeaks,
                                  function(x)metaData(x)$col))
groupNames <- unlist(lapply(mergedPeaks,
                            function(x)metaData(x)$group))
groupColors <- c("G_Obs"="#FF6600", "G_Pul"="#1B9AA5", "G_Nub"="#008040")
groupColors <- setNames(groupColors[groupNames], names(culicoidesColors))

colorLeafsSpecies <- function(x) {
  if (is.leaf(x)) {
    col <- culicoidesColors[attr(x, "label")]
    attr(x, "nodePar") <- list(lab.col=col, pch="")
    attr(x, "edgePar") <- list(col=col)
  }
  return(x)
}
colorLeafsGroup <- function(x) {
  if (is.leaf(x)) {
    col <- groupColors[attr(x, "label")]
    attr(x, "nodePar") <- list(lab.col=col, pch="")
    attr(x, "edgePar") <- list(col=col)
  }
  return(x)
}

msDendrogram <- dendrapply(as.dendrogram(msHclust), colorLeafsSpecies)
msGroupDendrogram <- dendrapply(msDendrogram, colorLeafsGroup)

## rotate C_Dew_2 and C_Dew_5 to the right extrema
msDendrogram[[2]][[2]][[2]] <- rev(msDendrogram[[2]][[2]][[2]])
msGroupDendrogram[[2]][[2]][[2]] <- rev(msGroupDendrogram[[2]][[2]][[2]])
## undo rotation for subtree
msDendrogram[[2]][[2]][[2]][[1]] <- rev(msDendrogram[[2]][[2]][[2]][[1]])
msGroupDendrogram[[2]][[2]][[2]][[1]] <- rev(msGroupDendrogram[[2]][[2]][[2]][[1]])

## plot dendrogram
postscript(file.path(plotPath, "ms_dendrogram.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
oldMar <- par("mar")
par(mar=c(2, 1, 1, 4), cex=0.5)
plot(msDendrogram, main="ms dendrogram", horiz=TRUE, center=FALSE)
bootlabels.hclust(as.hclust(msDendrogram), msBoot, horiz=TRUE, adj=c(-0.2, 0.5),
                  col=redgreen(msBoot))
par(mar=oldMar, cex=1)
dev.off()

postscript(file.path(plotPath, "heatmap.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
plotHeatmap(similarityMatrix, msDendrogram, msGroupDendrogram,
            speciesLabels=speciesLabels, groupLabels=groupLabels,
            speciesCol=culicoidesColors, groupCol=groupColors)
dev.off()

## plot pcr tree
pcrDendrogram <- dendrapply(as.dendrogram(pcrHclust), colorLeafsSpecies)
postscript(file.path(plotPath, "pcr_dendrogram.eps"),
           width=oneColumnWidth, height=oneColumnWidth)
oldMar <- par("mar")
par(mar=c(2, 1, 1, 7), cex=0.5)
plot(pcrDendrogram, main="pcr dendrogram", horiz=TRUE, center=FALSE, dLeaf=-0.012)
bootlabels.hclust(pcrHclust, pcrBoot, horiz=TRUE, adj=c(-0.2, 0.5),
                  col=redgreen(pcrBoot))
par(mar=oldMar, cex=1)
dev.off()

## plot pca
pch <- c(C_Nub=4, C_Dew=2, C_Obs=0, C_Sco=1,
         C_New=15, C_Pul=16, C_Pun=17)
pch <- setNames(pch[speciesLabels], names(speciesLabels))
postscript(file.path(plotPath, "pca.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plotPca(binaryPeakMatrix, speciesLabels, groupNames, culicoidesColors,
        groupColors, pch, conf=0.95)
dev.off()

## plot sda
postscript(file.path(plotPath, "sda_species.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plot(speciesRanking)
title(ylab="m/z")
dev.off()

postscript(file.path(plotPath, "sda_group.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plot(groupRanking)
title(ylab="m/z")
dev.off()

## plot half decimal place rule
postscript(file.path(plotPath, "hdpr.eps"),
           width=twoColumnWidth, height=oneColumnWidth)
plotHdpr(uniqueMass, hdprCo, tolerance=0.2)
dev.off()

###############################################################################
## END OF FIGURES
###############################################################################

