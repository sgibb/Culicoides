## Copyright 2012-2013 Sebastian Gibb
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

#' Averagine reference table
#'
#' This dataset contains mass for combinations of averagine.
#'
#' @details
#'  To create:
#'  \preformatted{
#'
#'averagineTable <- .createAveragineReferenceTable(averagine, minMass=120, maxMass=5000)
#'  }
#'
#' @references
#'  Michael W. Senko and Steven C. Beu and Fred W. McLafferty (1995), \cr
#'  \dQuote{Determination of monoisotopic masses and ion populations for large
#'  biomolecules from resolved isotopic distributions} \cr
#'  \emph{Journal of the American Society for Mass Spectrometry}, \bold{4} (6): 229-233 \cr
#'  ISSN: 1044-0305; doi 10.1016/1044-0305(95)00017-8 \cr
#'  \url{http://www.sciencedirect.com/science/article/pii/1044030595000178}
#'
#' @name averagineTable
#' @docType data
#' @rdname averagineTable
#' @keywords datasets

averagine <- data.frame(name="Averagine",
                        mass=111.1254,
                        C=4.9384,
                        H=7.7583,
                        N=1.3577,
                        O=1.4773,
                        S=0.0417, stringsAsFactors=FALSE)

.createAveragineReferenceTable <- function(averagine, minMass=120, maxMass=5000) {
  stopifnot(require("BRAIN"))

  n <- floor(minMass/averagine$mass) : ceiling(maxMass/averagine$mass)

  df <- data.frame()

  for (i in n) {
    l <- lapply(averagine[, c("C", "H", "N", "O", "S")], function(x)x*i)

    r <- .useBRAIN(l)

    if (!is.na(r$monoisotopicMass)) {
      df <- rbind(df, r)
    }
  }
  return(df)
}

#' This function is a wrapper around \code{\link[BRAIN]{useBRAIN}}.
#'
#' @param aC list with fields C, H, N, O, S of integer non-negative values.
#' @param nRatio number of ratio products.
#' @return data.frame
#'
#' @seealso \code{\link[BRAIN]{useBRAIN}}
#' @keywords internal
#' @rdname useBRAIN-functions
#'
.useBRAIN <- function(aC, nRatio=2) {
  r <- BRAIN::useBRAIN(aC=aC, nrPeaks=1000, stopOption="abundantEstim",
                       abundantEstim=10)

  maxIdx <- which.max(r$isoDistr)
  if (is.na(r$monoisotopicMass)) {
    return(data.frame())
  } else {
    d <- data.frame(monoisotopicMass=r$monoisotopicMass,
                    apexMass=r$masses[maxIdx],
                    apexIdx=maxIdx,
                    MvsA=r$isoDistr[1]/r$isoDistr[maxIdx])
    s <- 1:(nRatio+1)
    r$isoDistr <- r$isoDistr[s]
    q <- r$isoDistr[-1]/head(r$isoDistr, -1)
    d <- cbind(d, t(q))
    s <- s-1
    colnames(d) <- c("monoisotopicMass", "apexMass", "apexIdx", "MvsA",
                     paste("I", s[-1], "vsI", head(s, -1), sep=""))
    return(as.data.frame(d))
  }
}

