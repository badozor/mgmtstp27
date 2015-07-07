## Title ##
Prediction of DNA methylation state of MGMT promoter based on HM-450K and HM-27K infinium platforms

## Description ##
This R package contains function to compute the prediction of the DNA methylation of MGMT promoter with data from infinium HM-450K and HM-27K platforms

## License ##
GPL version 2 or newer
```
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
```

## Installation ##
To install this package, load the archive containing the package mgmtstp27 (see below), start R and enter:

```
install.packages("mgmtstp27_0.6.tar.gz",repos=NULL)
install.packages("mgmtstp27_0.6.zip",repos=NULL)
```

The sources are avalaible [here](https://github.com/badozor/mgmtstp27/tree/master/trunk/Rpackage). 

## Example ##
The R function MGMTpredict directly provides prediction, classification and confidence intervals as illustrated below:
```
# loading R packages
require(mgmtstp27)
require(minfiData)
# preprocessing of the data
dat <- preprocessRaw(RGsetEx)
# computation of M-value
mvalue <- log2((getMeth(dat)+1)/(getUnmeth(dat)+1))
mvalue <- as.data.frame(t(mvalue))
# predictions
pred1 <- MGMTpredict(mvalue)
head(pred1)
# quality control graphics
par(mfrow=c(2,3))
MGMTqc.pop(pred1,which.plot=1:3,mfrow=NULL)
MGMTqc.single(pred1,nsample=1,which.plot=1:3,mfrow=NULL)
```
Additional documentations are avaialble [here](https://github.com/badozor/mgmtstp27/tree/master/trunk/Rdoc):
  * Introduction to R package mgmtstp27 (version 0.1, in preparation, intromgmtstp27.pdf)
  * Effect of normalization on the prediction of DNA methylation status of MGMT promoter: example with HM-450K Infinium data from TCGA and the R package mgmtstp27 (version 0.1, in preparation, docmgmtstp27.pdf)
  * Prediction of the DNA methylation of MGMT with raw data (format IDAT) from HM-27k (version 0.1, in preparation, MgmtPredTCGA.pdf)


## Unexhaustive References ##
  * Bady, P., D. Sciuscio, A.-C. Diserens, J. Bloch, M. J. van den Bent, C. Marosi, P.-Y. Dietrich, M. Weller, L. Mariani, F. L. Heppner, D. R. McDonald, D. Lacombe, R. Stupp, M. Delorenzi, and M. E. Hegi. (2012). MGMT methylation analysis of glioblastoma on the Infinium methylation BeadChip identifies two distinct CpG regions associated with gene silencing and outcome, yielding a prediction model for comparisons across datasets, tumor grades, and CIMP-status. Acta Neuropathologica 124:547-560. PubMed:http://www.ncbi.nlm.nih.gov/pubmed/22810491

## Depends ##
R (>= 3.1.2), minfi, lumi, ade4,methylumi,MASS


R (>= 3.2.0), minfi, lumi, ade4,methylumi,MASS

## Suggests ##
boot

## Date ##
2014-09-11

## Revison ##
2015-06-16

## Version ##
alpha 0.6 (version for R-3.1.2)


alpha 0.6-2 (version for R-3.2.0)

## URL ##
  * http://lausanne.isb-sib.ch/~pbady/Rpackages.html
  * http://www.chuv.ch/neurosciences/en/dnc-recherche-laboratoire_de_biologie_et_genetique_des_tumeurs_cerebrales.htm