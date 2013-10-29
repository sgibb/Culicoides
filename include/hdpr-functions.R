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

decimalPlace <- function(x) {
  return(x-trunc(x))
}

decimalPlaceHdpr <- function(x) {
  d <- decimalPlace(x)
  return(d + round((x * 5e-4 - d)))
}

nominalMass <- function(x) {
  return(x - decimalPlaceHdpr(x))
}

## Calculate HDPR regression for mass x
## Needed for filtering using HDPR (Half-Decimal-Place-Rule)
hdprCoef <- function(x) {
  rlm <- MASS::rlm(y ~ x, data=list(y=decimalPlaceHdpr(x), x=nominalMass(x)))
  return(coef(rlm))
}

## filter by HDPR
## params:
##  x: double, mass values
##  co: coefficients; each element must be a double vector with two
##  elements (output of rlm)
##  tolerance: +/-tolerance
## returns:
##  logical (length == length(x))
##  TRUE: in HDPR range, FALSE: out of HDPR range
hdprFilter <- function(x, co, tolerance=0.2) {
  d <- decimalPlaceHdpr(x)
  m <- nominalMass(x)

  upper <- (co[1]+tolerance) + co[2]*m
  lower <- (co[1]-tolerance) + co[2]*m

  return(lower <= d & d <= upper)
}

