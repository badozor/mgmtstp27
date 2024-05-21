## Title ##
Prediction of DNA methylation state of MGMT promoter based on infinium platforms (27K, 450K, EPIC and EPICv2)

## Description ##
This R package contains function to compute the prediction of the DNA methylation of MGMT promoter with data from infinium EPIC, EPIC version 2, HM-450K and HM-27K platforms

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
install.packages("mgmtstp27_x.x-x.tar.gz",repos=NULL)
install.packages("mgmtstp27_x.x-x.zip",repos=NULL)
```

where x.x-x corresponds to the version of the package.

The sources are avalaible [here](https://github.com/badozor/mgmtstp27/tree/master/trunk/Rpackage). 


The installation of the package mgmtstp27 requires the presence (or installation) of the following additional packages:
```
## CRAN (updated 21/05/2024)
install.packages(c("ade4","MASS"))

## Bioconductor (https://www.bioconductor.org/install/)
## try http:// if https:// URLs are not supported
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("lumi","methylumi","minfi"))
```


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


## Note for EPIC v2 infinium platform (update 18/05/2024)

For EPIC v2 infinium platform, the function uses the two probes "cg12434587_BC11" and "cg12981137_TC11" (which will be renamed as "cg12434587" and "cg12981137"). The selection of these two probes is based on the concordance of the annotation information and the sequences of the probes (see variable "AlleleA_ProbeSeq" from manifest file). The explanation is given below:


```
# import R packages
require(mgmtstp27)
coef(MGMTSTP27)
probemodels <- names(coef(MGMTSTP27))[-1]

# data source: https://zwdzwd.github.io/InfiniumAnnotation
# manifest for EPIC v2
annotepicv2 <- read.table("EPICv2.hg38.manifest.tsv",header=TRUE,sep="\t")
dim(annotepicv2)
head(annotepicv2)
probenamesv2 <- annotepicv2$Probe_ID

# manifest for EPIC
annotepic <- read.table("EPIC.hg38.manifest.tsv",header=TRUE,sep="\t")
dim(annotepic)
head(annotepic)
probenames <- annotepic$Probe_ID

##
# checknames
V2probes1 <- annotepicv2[grep(probemodels[1],probenamesv2),]
V2probes2 <- annotepicv2[grep(probemodels[2],probenamesv2),]
V1probes1 <- annotepic[grep(probemodels[1],probenames),]
V1probes2 <- annotepic[grep(probemodels[2],probenames),]

## for probe 1: "cg12434587"

V1probes1
V2probes1
V2probes1[is.element(V2probes1$AlleleA_ProbeSeq,V1probes1$AlleleA_ProbeSeq),]

## for probe 2: "cg12981137"
V1probes2
V2probes2
V2probes2[is.element(V2probes2$AlleleA_ProbeSeq,V1probes2$AlleleA_ProbeSeq),]

## final selection
V2probes1[is.element(V2probes1$AlleleA_ProbeSeq,V1probes1$AlleleA_ProbeSeq),"Probe_ID"]
V2probes2[is.element(V2probes2$AlleleA_ProbeSeq,V1probes2$AlleleA_ProbeSeq),"Probe_ID"]
```


## Depends ##
R (>= 3.1.2), minfi, lumi, ade4,methylumi,MASS

R (>= 3.2.0), minfi, lumi, ade4,methylumi,MASS

R (>= 3.2.2), minfi, lumi, ade4,methylumi,MASS

R (>= 4.2.1), minfi, lumi, ade4,methylumi,MASS

## Suggests ##
boot

## Date ##
2014-09-11

## Revison ##
2022-07-22

## Version ##
0.6 (version for R-3.1.2)

0.6-2 (version for R-3.2.0)

0.6-3 (version for R-3.2.2)

0.7 (version for R-4.2.1)

0.8 (version for R-4.4.0)

## URL ##
  * http://lausanne.isb-sib.ch/~pbady/Rpackages.html
  * http://www.chuv.ch/neurosciences/en/dnc-recherche-laboratoire_de_biologie_et_genetique_des_tumeurs_cerebrales.htm


## Unexhaustive References ##
  * Bady, P., D. Sciuscio, A.-C. Diserens, J. Bloch, M. J. van den Bent, C. Marosi, P.-Y. Dietrich, M. Weller, L. Mariani, F. L. Heppner, D. R. McDonald, D. Lacombe, R. Stupp, M. Delorenzi, and M. E. Hegi. (2012). MGMT methylation analysis of glioblastoma on the Infinium methylation BeadChip identifies two distinct CpG regions associated with gene silencing and outcome, yielding a prediction model for comparisons across datasets, tumor grades, and CIMP-status. Acta Neuropathologica 124:547-560. PubMed:http://www.ncbi.nlm.nih.gov/pubmed/22810491
  * van den Bent MJ, Erdem-Eraslan L, Idbaih A, de Rooi J, Eilers PHC, Spliet WGM, den Dunnen WFA, Tijssen C, Wesseling P, Sillevis Smitt PAEet al (2013) MGMT-STP27 Methylation Status as Predictive Marker for Response to PCV in Anaplastic Oligodendrogliomas and Oligoastrocytomas. A Report from EORTC Study 26951. Clinical Cancer Research 19: 5513-5522. PubMed: http://www.ncbi.nlm.nih.gov/pubmed/23948976
  * Bady P., Delorenzi M., Hegi M. (accepted) Sensitivity analysis of the MGMT-STP27 model and impact of genetic/epigenetic context to predict the MGMT methylation status in gliomas and other tumors, Journal of Molecular Diagnostics, xx,xxxx-xxxx.PubMed: http://www.ncbi.nlm.nih.gov/pubmed/26927331
  * Hench IB, Monica RD, Chiariotti L, Bihl M, Tolnay M, Frank S, Hench J (2021) Fast routine assessment of MGMT promoter methylation. Neurooncol Adv 3: vdaa170 Doi 10.1093/noajnl/vdaa170
  * Tzaridis T, Schafer N, Weller J, Steinbach JP, Schlegel U, Seidel S, Sabel M, Hau P, Seidel C, Krex Det al (2021) MGMT promoter methylation analysis for allocating combined CCNU/TMZ chemotherapy: Lessons learned from the CeTeG/NOA-09 trial. Int J Cancer 148: 1695-1707 Doi 10.1002/ijc.3336
  * zzz
 

## A list of paper citing the model MGMT-STP27 ##

Cited In for PMID: 22810491 from pubmed: [here](https://pubmed.ncbi.nlm.nih.gov/?size=200&linkname=pubmed_pubmed_citedin&from_uid=22810491)





