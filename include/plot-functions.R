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

plotSpectra <- function(spectra, main, xlim, ylim,
                        speciesMass, groupMass,
                        yaxt="n", lwd=0.25,
                        speciesPeaks.pch=8, groupPeaks.pch=6,
                        peaks.lwd=lwd, peaks.cex=0.5, add.expr) {
  n <- length(spectra)

  oldPar <- par(no.readonly=TRUE)
  par(mfrow=c(n, 1))

  if (missing(main)) {
    main <- unlist(lapply(spectra, function(x)metaData(x)$name))
  }

  if (missing(xlim)) {
    xlim <- range(unlist(lapply(spectra, function(x)range(mass(x)))))
  }

  par(oma=c(2, 1, 0, 1), mar=rep(0, 4), bty="o", xaxs="i", tcl=-0.2,
      mgp=c(3, 0.25, 0))

  for (i in seq(along=spectra)) {
    m <- main[i]
    col <- metaData(spectra[[i]])$col
    if (is.null(col)) {
      col <- 1
    }

    if (missing(ylim)) {
      massRange <- mass(spectra[[i]]) >= xlim[1] &
                   mass(spectra[[i]]) <= xlim[2]
      mI <- max(intensity(spectra[[i]])[massRange])
      yli <- c(0, mI*1.01)
    } else {
      if (all(ylim >= 0 & ylim <= 2)) {
        yli <- ylim*max(spectra[[i]]@intensity)
      } else {
        yli <- ylim
      }
    }

    plot(spectra[[i]], xlim=xlim, ylim=yli,
         main="", sub="", xaxt="n", yaxt=yaxt, ylab="",
         las=1, frame.plot=FALSE, type="n")
    if (!missing(speciesMass)) {
      abline(v=unique(as.vector(speciesMass)), col="#808080", lwd=peaks.lwd*2,
             lty=4)
    }
    if (!missing(groupMass)) {
      abline(v=unique(as.vector(groupMass)), col="#808080", lwd=peaks.lwd,
             lty=5)
    }
    lines(spectra[[i]], col=col, lwd=lwd)
    box(lwd=0.25)

    if (i == n) {
      labels <- tail(head(pretty(xlim), -1), -1)
      at <- c(xlim[1], labels, xlim[2])
      axis(side=1, at=at, labels=c(character(1), labels, character(1)))
    }

    if (!missing(speciesMass)) {
      points(speciesMass[, i], rep(yli[2]*0.99, nrow(speciesMass)), col=col,
             pch=speciesPeaks.pch, lwd=peaks.lwd, cex=peaks.cex)
    }
    if (!missing(groupMass)) {
      points(groupMass[, metaData(spectra[[i]])$group],
             rep(yli[2]*0.99, nrow(groupMass)), col=col,
             pch=groupPeaks.pch, lwd=peaks.lwd, cex=peaks.cex)
    }

    maxStrWidth <- max(strwidth(main), na.rm=TRUE)
    maxStrHeight <- max(strheight(main), na.rm=TRUE)

    usr <- par("usr")
    x2 <- usr[2]+diff(usr[1:2])/1.08*(-0.02)-maxStrWidth
    y2 <- usr[4]+diff(usr[3:4])/1.08*(-0.04)
    text(x2, y2, labels=m, adj=c(0, 1), col=1, font=3) # italic
    text(x2, y2-maxStrHeight*1.5, labels=metaData(spectra[[i]])$spot,
         adj=c(0, 1), col="#808080", font=1, cex=0.5)
  }
  title(xlab="m/z", outer=TRUE, line=1)
  title(ylab="normalized intensity", outer=TRUE, line=0)

  if (!missing(add.expr)) {
    eval(substitute(add.expr))
  }

  par(oldPar)
}

plotZoomedSpectra <- function(spectra, main, pos, delta, xlim, ylim,
                              speciesMass, groupMass,
                              yaxt="n", lwd=0.25,
                              speciesPeaks.pch=8, groupPeaks.pch=6,
                              peaks.lwd=lwd, peaks.cex=0.5) {
  ns <- length(spectra)
  np <- length(pos)
  sortedMass <- sort(c(speciesMass, groupMass))

  oldPar <- par(no.readonly=TRUE)
  par(mfcol=c(ns, np))

  if (missing(main)) {
    main <- unlist(lapply(spectra, function(x)metaData(x)$name))
  }

  par(oma=c(2, 1, 0, 1), mar=c(rep(0, 3), 0.25), bty="o", xaxs="i", tcl=-0.2,
      mgp=c(3, 0.25, 0))

  for (j in seq(along=pos)) {
    for (i in seq(along=spectra)) {
      m <- main[i]
      col <- metaData(spectra[[i]])$col
      if (is.null(col)) {
        col <- 1
      }

      xlim <- pos[j]+c(-delta, +delta)

      if (missing(ylim)) {
        massRange <- mass(spectra[[i]]) >= xlim[1] &
                     mass(spectra[[i]]) <= xlim[2]
        mI <- max(intensity(spectra[[i]])[massRange])
        yli <- c(0, mI+0.01)
      } else {
        if (all(ylim >= 0 & ylim <= 2)) {
          yli <- ylim*max(spectra[[i]]@intensity)
        } else {
          yli <- ylim
        }
      }

      plot(spectra[[i]], xlim=xlim, ylim=yli,
           main="", sub="", xaxt="n", yaxt=yaxt, ylab="",
           las=1, frame.plot=FALSE, type="n")
      if (!missing(speciesMass)) {
        abline(v=unique(as.vector(speciesMass)), col="#808080", lwd=peaks.lwd,
               lty=4)
      }
      if (!missing(groupMass)) {
        abline(v=unique(as.vector(groupMass)), col="#808080", lwd=peaks.lwd,
               lty=5)
      }
      lines(spectra[[i]], col=col, lwd=lwd)
      box(lwd=0.25)

      if (i == ns) {
        labels <- round(sortedMass[MALDIquant:::.which.closest(pos[j],
                                                               sortedMass)])
        at <- c(xlim[1], labels, xlim[2])
        axis(side=1, at=at, labels=c(character(1), labels, character(1)))
      }

      if (!missing(speciesMass)) {
        points(speciesMass[, i], rep(yli[2]*0.99, nrow(speciesMass)), col=col,
               pch=speciesPeaks.pch, lwd=peaks.lwd, cex=peaks.cex)
      }
      if (!missing(groupMass)) {
        points(groupMass[, metaData(spectra[[i]])$group],
               rep(yli[2]*0.99, nrow(groupMass)), col=col,
               pch=groupPeaks.pch, lwd=peaks.lwd, cex=peaks.cex)
      }
    }
  }
  title(xlab="m/z", outer=TRUE, line=1)
  title(ylab="normalized intensity", outer=TRUE, line=0)
  par(oldPar)
}

plotHeatmap <- function(x, dendrogram, groupDendrogram,
                        speciesLabels, groupLabels,
                        speciesCol, groupCol) {
  col <- colorRampPalette(c("#000000", "#FF0000", "#FFFF00"))(100)

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  heatmap(x, Colv=dendrogram, Rowv=groupDendrogram,
          symm=TRUE, revC=TRUE, margins=c(10, 10),
          ColSideColor=speciesCol[colnames(x)],
          RowSideColor=groupCol[colnames(x)],
          col=col,
          add.expr=abline(v=seq(0.5, ncol(x)+0.5, length.out=ncol(x)+1),
                          h=seq(0.5, ncol(x)+0.5, length.out=ncol(x)+1),
                          col="#808080", lwd=0.25))

  text(x=0.95, y=0.95, "species", cex=0.5, font=2)
  legend(x=0.91, y=0.94, legend=unique(speciesLabels),
         col=unique(speciesCol), pch=15, cex=0.5, bty="n")

  text(x=0.95, y=0.70, "groups", cex=0.5, font=2)
  legend(x=0.91, y=0.69, legend=unique(groupLabels),
         col=unique(groupCol), pch=15, cex=0.5, bty="n")

  #par(new=TRUE, pty="m", plt=c(0.92, 0.94, 0.1, 0.25), cex=0.5, err=-1,
  par(new=TRUE, pty="m", plt=c(0.92, 0.94, 0.25, 0.4), cex=0.5, err=-1,
      mgp=c(0, 0.2, 0.2))
  image(matrix(1:100, nrow=1), col=col, xaxt="n", yaxt="n")
  axis(2, at=seq(0, 1, by=0.5), tick=FALSE, las=1, line=0)
}

plotPhyloTree <- function(tree, bootstrap, col, type="phylogram", threshold=80,
                          cex=0.6, ...) {
  mar <- par("mar")
  on.exit(par(mar=mar, cex=1))
  par(mar=c(4, 0, 0, 0), cex=cex)
  ape::plot.phylo(tree, type=type, cex=cex*1.5, tip.color=col[tree$tip.label],
                  label.offset=0.0075, underscore=TRUE)
  axis(1)
  b <- round(bootstrap/max(bootstrap)*100, 0)
  par(cex=cex*0.75)
  nodelabels(b, col=ifelse(b >= threshold, "#4DAF4A", "#E41A1C"),
             frame="none", adj=-0.1)
}

plotPca <- function(x, speciesNames, groupNames, individualCol, groupCol,
                    pch, conf=0.95) {
  usN <- unique(speciesNames)
  ugN <- unique(groupNames)
  ugC <- setNames(unique(groupCol), ugN)

  p <- vector(mode="list", length=length(ugN)+1)

  pca <- .pca(x, groups=groupNames, species=substr(rownames(x), 1, 5),
              conf=conf)
  p[[1]] <- list(pca=pca, pch=pch, pcol=individualCol,
                 ecol=tapply(groupCol, groupNames,
                             function(x)x[1])[names(pca$e)])

  p[-1] <- tapply(rownames(x), groupNames, function(y) {
    if (length(unique(speciesNames[y])) > 1) {
      pca <- .pca(x[y, ], groups=speciesNames[y], conf=conf)
      return(list(pca=pca,
                  pch=pch[y], pcol=individualCol[y],
                  ecol=tapply(individualCol[y], speciesNames[y],
                              function(x)x[1])[names(pca$e)]))
    } else {
      return(NULL)
    }
  })

  oldPar <- par(no.readonly=TRUE)
  par(mar=rep(0, 4), tcl=-0.2)
  isValidPca <- !unlist(lapply(p, is.null))
  nr <- sum(isValidPca)
  plotLabels <- character(nr)
  plotLabels[isValidPca] <- LETTERS[1:nr]
  m <- matrix(c(rep(1, nr-1), 2:nr), byrow=FALSE, ncol=2)
  layout(m, widths=c(2, 1))
  for (i in seq(along=p)) {
    if (!is.null(p[[i]])) {
      .plotPca(p[[i]]$pca$p, hull=p[[i]]$pca$h, ellipse=p[[i]]$pca$e,
               ellipse.label=(i==1), plot.label=plotLabels[i],
               pch=p[[i]]$pch, pcol=p[[i]]$pcol, ecol=p[[i]]$ecol)
      if (i == 1) {
        ord <- unlist(lapply(split(speciesNames, groupNames),
                             function(x)sort(unique(x))))
        sel <- which(!duplicated(speciesNames))
        sel <- sel[match(ord, speciesNames[sel])]
        legend("topleft", legend=speciesNames[sel],
               col=individualCol[sel], pch=pch[sel], bty="n")
      }
    }
  }
  par(oldPar)
  invisible(p)
}

.pca <- function(x, groups, species=groups, conf) {
  p <- vegan::rda(x, scale=FALSE)
  h <- vegan::ordihull(p, display="sites", draw="none",
                       groups=species)
  e <- vegan::ordiellipse(p, display="sites", draw="none", groups=groups,
                          conf=conf)
  return(list(p=p, h=h, e=e))
}

.plotPca <- function(p, hull, ellipse,
                     ellipse.label, plot.label,
                     pch, pcol, ecol,
                     lty.hull=1, lty.ellipse=2,
                     lwd.hull=1, lwd.ellipse=1) {
  lab <- paste("Component ", 1:2, " (", round(p$CA$eig[1:2], digits=2),
               " %)", sep="")
  par(mar=c(3.5, 3, 0, 1))
  plot(p, type="n", display="sites", xlab="", ylab="")
  .labelPlot(plot.label)
  title(xlab=lab[1], line=2)
  title(ylab=lab[2], line=2)

  points(p, display="sites", pch=pch, col=pcol)

  for (i in seq(along=hull)) {
    lines(hull[[i]], lwd=lwd.hull, lty=lty.hull,
          col=pcol[rownames(hull[[i]])[1]])
  }

  needLabel <- names(ellipse) != "G_Nub"
  ellipseLabels <- character(length(ellipse))
  ellipseLabels[needLabel] <- letters[(1:sum(needLabel))+1]

  for (i in seq(along=ellipse)) {
    xy <- vegan:::veganCovEllipse(ellipse[[i]]$cov, ellipse[[i]]$center,
                                  ellipse[[i]]$scale)
    if (ellipse.label) {
      text(ellipse[[i]]$center[1], ellipse[[i]]$center[2],
           labels=ellipseLabels[i], cex=1.25, font=3)
    }
    lines(xy, lwd=lwd.ellipse, lty=lty.ellipse, col=ecol[i])
  }
}

.labelPlot <- function(char, cex=1.5) {
    usr <- par("usr")
    text(x=usr[2]-(cex*strwidth(char)),
         y=usr[4]-(cex*strheight(char)),
         font=2, labels=char, cex=cex)
}

plotHdpr <- function(x, co=c(0, 5e-4), tolerance=0.2) {
  ## plot peak mass vs decimal place
  decimalPlaces <- decimalPlaceHdpr(x)

  f <- hdprFilter(x, co=co, tolerance=tolerance)
  m <- nominalMass(x)

  plot(m[f], decimalPlaces[f], main="HDPR",
       xlab="mass", ylab="decimal place",
       ylim=c(0, ceiling(max(decimalPlaces))))
  mtext(paste0("y=", round(co[2], 5), "*x+", round(co[1], 5)))

  points(m[!f], decimalPlaces[!f], col=1, pch=4)

  abline(coef=c(co[1]-tolerance, co[2]), col=2, lty=2)
  abline(coef=c(co[1], co[2]), col=2)
  abline(coef=c(co[1]+tolerance, co[2]), col=2, lty=2)

  legend("bottomright", legend=c("monoisotopic mass", "filtered mass",
                                 "HDPR", "HDPR tolerance range"),
         col=rep(1:2, each=2), pch=c(1, 4, NA, NA),
         lwd=1, lty=c(NA, NA, 1, 2), bg="white")
}

redgreen <- function(x, threshold=0.8) {
  return(ifelse(x < threshold, "red", "darkgreen"))
}

