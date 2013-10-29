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

topSdaPeaks <- function(speciesRanking, groupRanking, xlim,
                        top=c(species=3, group=3)) {
  if (missing(xlim)) {
    xlim <- range(as.double(rownames(speciesRanking)))
  }

  .topPeaks <- function(x, xlim, top, exclude=NA) {
    m <- as.double(rownames(x))
    cl.count <- attr(x, "cl.count")
    isDDA <- attr(x, "diagonal")
    idx <- xlim[1] <= m & m <= xlim[2] & !m %in% exclude
    m <- m[idx]
    x <- x[idx, (1:cl.count)+2]
    x <- apply(x, 2, function(y) {
      ix <- sort(y, method="quick", decreasing=TRUE,
                 index.return=TRUE)$ix[1:top]
      return(as.double(m[ix]))
    })
    colnames(x) <- substring(colnames(x), ifelse(isDDA, 3, 5))
    return(x)
  }

  groupRanking <- .topPeaks(groupRanking, xlim, top["group"])

  speciesRanking <- .topPeaks(speciesRanking, xlim, top["species"],
                              groupRanking)

  return(list(species=speciesRanking, group=groupRanking))
}
