## Title ##
Prediction of DNA methylation state of MGMT promoter based on HM-450K and HM-27K infinium platforms

## Description ##
This R package contains function to compute the prediction of the DNA methylation of MGMT promoter with data from infinium HM-450K and HM-27K platforms


**NOTA:** This version is a beta version! as results, the R package mgmtstp27 is still in development and the documentation need some corrections (e.g. language, etc ...).


## Installation ##
To install this package, load the archive containing the package mgmtstp27 (see below), start R and enter:

```
install.packages("mgmtstp27_0.6.tar.gz",repos=NULL)
install.packages("mgmtstp27_0.6.zip",repos=NULL)
```


The sources are avalaible [here](https://code.google.com/p/mgmtstp27/source/browse/#git%2Ftrunk%2FRpackage).


Additional documentation is avaialble [here](https://code.google.com/p/mgmtstp27/source/browse/#git%2Ftrunk%2FRdoc):
  * Introduction to R package mgmtstp27 (version 0.1, in preparation, intromgmtstp27.pdf)
  * Effect of normalization on the prediction of DNA methylation status of MGMT promoter: example with HM-450K Infinium data from TCGA and the R package mgmtstp27 (version 0.1, in preparation, docmgmtstp27.pdf)
  * Prediction of the DNA methylation of MGMT with raw data (format IDAT) from HM-27k (version 0.1, in preparation, MgmtPredTCGA.pdf)


## Unexhaustive References ##
  * Bady, P., D. Sciuscio, A.-C. Diserens, J. Bloch, M. J. van den Bent, C. Marosi, P.-Y. Dietrich, M. Weller, L. Mariani, F. L. Heppner, D. R. McDonald, D. Lacombe, R. Stupp, M. Delorenzi, and M. E. Hegi. (2012). MGMT methylation analysis of glioblastoma on the Infinium methylation BeadChip identifies two distinct CpG regions associated with gene silencing and outcome, yielding a prediction model for comparisons across datasets, tumor grades, and CIMP-status. Acta Neuropathologica 124:547-560. PubMed:http://www.ncbi.nlm.nih.gov/pubmed/22810491

## Depends ##
R (>= 3.1.2), minfi, lumi, ade4,methylumi,MASS

## Suggests ##
boot

## License ##
GPL version 2 or newer

## Date ##
2014-09-11

## Revison ##
2015-06-05

## Version ##
alpha 0.6 (in construction)

## URL ##
  * http://lausanne.isb-sib.ch/~pbady/Rpackages.html
  * http://www.chuv.ch/neurosciences/en/dnc-recherche-laboratoire_de_biologie_et_genetique_des_tumeurs_cerebrales.htm