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

## find representative spots (individual spectra) by look for the highest mean
## similarity score
representativeSpots <- function(peakList, speciesLabels) {
  n <- length(peakList)

  if (missing(speciesLabels)) {
    speciesLabels <- unlist(lapply(peakList, function(x) {
      return(metaData(x)$species)
    }))
  }

  binaryMatrix <- ifelse(is.na(intensityMatrix(peakList)), 0, 1)
  similarityMatrix <- as.matrix(proxy::simil(binaryMatrix, method="Dice"))

  s <- tapply(1:n, speciesLabels, function(i) {
    m <- rowMeans(similarityMatrix[i, i], na.rm=TRUE)
    return(i[which.max(m)])
  })

  return(s)
}

