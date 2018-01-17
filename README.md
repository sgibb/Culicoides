# Species determination of Culicoides biting midges via peptide profiling using matrix-assisted laser desorption ionization mass spectrometry

## Introduction

This repository contains the complete [R](http://r-project.org) scripts to reproduce the analysis described in:

K. R. Uhlmann, S.Gibb, S. Kalkhof, U. Arroyo-Abad, C. Schulz,
B. Hoffmann, F. Stubbins, S. Carpenter, M. Beer, M. von Bergen, R. Feltens. 2013.
*Species determination of Culicoides biting midges via peptide profiling using
matrix-assisted laser desorption ionization mass spectrometry*. https://doi.org/10.1186/1756-3305-7-392

Please find the dataset at http://dx.doi.org/10.6084/m9.figshare.801878 <br />
(The dataset is downloaded by the `setup.R` script automatically.)

## Install

 - You have to install R (>= 3.0.0) for your operation system: http://cran.r-project.org/
 - Download the R scripts: https://github.com/sgibb/Culicoides/archive/master.zip
 - Start R and run the *setup.R* script (this is needed only for the first use):

```s
source("setup.R")
```

## Run the analysis

After installing R and the needed packages by running the *setup.R* script you could start the analysis by calling:

```s
source("analysis.R")
```

(use `source("figure7")` to create the *IPP vs SMM* plot (figure 7) and `source("digest_drosophila_melanogaster.R")` to perform the in-silico digestion of the *Drosophila melanogaster* reference proteome)

## Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/sgibb/Culicoides/issues>
* compose an e-mail to: <mail@sebastiangibb.de>

