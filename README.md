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
install.packages("mgmtstp27_x.x-x.tar.gz",repos=NULL)
install.packages("mgmtstp27_x.x-x.zip",repos=NULL)
```

where x.x-x corresponds to the version of the package.

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


## Depends ##
R (>= 3.1.2), minfi, lumi, ade4,methylumi,MASS

R (>= 3.2.0), minfi, lumi, ade4,methylumi,MASS

R (>= 3.2.2), minfi, lumi, ade4,methylumi,MASS

## Suggests ##
boot

## Date ##
2014-09-11

## Revison ##
2015-11-26

## Version ##
0.6 (version for R-3.1.2)

0.6-2 (version for R-3.2.0)

0.6-3 (version for R-3.2.2)

## URL ##
  * http://lausanne.isb-sib.ch/~pbady/Rpackages.html
  * http://www.chuv.ch/neurosciences/en/dnc-recherche-laboratoire_de_biologie_et_genetique_des_tumeurs_cerebrales.htm


## Unexhaustive References ##
  * Bady, P., D. Sciuscio, A.-C. Diserens, J. Bloch, M. J. van den Bent, C. Marosi, P.-Y. Dietrich, M. Weller, L. Mariani, F. L. Heppner, D. R. McDonald, D. Lacombe, R. Stupp, M. Delorenzi, and M. E. Hegi. (2012). MGMT methylation analysis of glioblastoma on the Infinium methylation BeadChip identifies two distinct CpG regions associated with gene silencing and outcome, yielding a prediction model for comparisons across datasets, tumor grades, and CIMP-status. Acta Neuropathologica 124:547-560. PubMed:http://www.ncbi.nlm.nih.gov/pubmed/22810491
  * van den Bent MJ, Erdem-Eraslan L, Idbaih A, de Rooi J, Eilers PHC, Spliet WGM, den Dunnen WFA, Tijssen C, Wesseling P, Sillevis Smitt PAEet al (2013) MGMT-STP27 Methylation Status as Predictive Marker for Response to PCV in Anaplastic Oligodendrogliomas and Oligoastrocytomas. A Report from EORTC Study 26951. Clinical Cancer Research 19: 5513-5522. PubMed: http://www.ncbi.nlm.nih.gov/pubmed/23948976
  * Bady P., Delorenzi M., Hegi M. (accepted) Sensitivity analysis of the MGMT-STP27 model and impact of genetic/epigenetic context to predict the MGMT methylation status in gliomas and other tumors, Journal of Molecular Diagnostics, xx,xxxx-xxxx.PubMed: http://www.ncbi.nlm.nih.gov/pubmed/26927331
  * zzz
 

## A list of paper citing the model MGMT-STP27 ##
last update: 2015-07-09 from pubmed

  * Wesseling P, van den Bent M, Perry A. Oligodendroglioma: pathology, molecular 
mechanisms and markers. Acta Neuropathol. 2015 Jun;129(6):809-27. doi:
10.1007/s00401-015-1424-1. Epub 2015 May 6. PubMed PMID: 25943885; PubMed Central
PMCID: PMC4436696.
  * Vigneswaran K, Neill S, Hadjipanayis CG. Beyond the World Health Organization 
grading of infiltrating gliomas: advances in the molecular genetics of glioma
classification. Ann Transl Med. 2015 May;3(7):95. doi:
10.3978/j.issn.2305-5839.2015.03.57. Review. PubMed PMID: 26015937; PubMed
Central PMCID: PMC4430738.

  * van Thuijl HF, Mazor T, Johnson BE, Fouse SD, Aihara K, Hong C, Malmström A,
Hallbeck M, Heimans JJ, Kloezeman JJ, Stenmark-Askmalm M, Lamfers ML, Saito N,
Aburatani H, Mukasa A, Berger MS, Söderkvist P, Taylor BS, Molinaro AM, Wesseling
P, Reijneveld JC, Chang SM, Ylstra B, Costello JF. Evolution of DNA repair
defects during malignant progression of low-grade gliomas after temozolomide
treatment. Acta Neuropathol. 2015 Apr;129(4):597-607. doi:
10.1007/s00401-015-1403-6. Epub 2015 Feb 28. PubMed PMID: 25724300; PubMed
Central PMCID: PMC4482618.

  * Kurscheid S, Bady P, Sciuscio D, Samarzija I, Shay T, Vassallo I, Criekinge
WV, Daniel RT, van den Bent MJ, Marosi C, Weller M, Mason WP, Domany E, Stupp R, 
Delorenzi M, Hegi ME. Chromosome 7 gain and DNA hypermethylation at the HOXA10
locus are associated with expression of a stem cell related HOX-signature in
glioblastoma. Genome Biol. 2015 Jan 27;16:16. doi: 10.1186/s13059-015-0583-7.
PubMed PMID: 25622821; PubMed Central PMCID: PMC4342872.

  * Le Rhun E, Taillibert S, Chamberlain MC. The future of high-grade glioma:
Where we are and where are we going. Surg Neurol Int. 2015 Feb 13;6(Suppl
1):S9-S44. doi: 10.4103/2152-7806.151331. eCollection 2015. Erratum in: Surg
Neurol Int. 2015 Mar 5;6:37. Rhun, Emilie Le [corrected to Le Rhun, Emilie].
PubMed PMID: 25722939; PubMed Central PMCID: PMC4338495.

  * Wiestler B, Capper D, Hovestadt V, Sill M, Jones DT, Hartmann C, Felsberg J,
Platten M, Feiden W, Keyvani K, Pfister SM, Wiestler OD, Meyermann R,
Reifenberger G, Pietsch T, von Deimling A, Weller M, Wick W. Assessing CpG island
methylator phenotype, 1p/19q codeletion, and MGMT promoter methylation from
epigenome-wide data in the biomarker cohort of the NOA-04 trial. Neuro Oncol.
2014 Dec;16(12):1630-8. doi: 10.1093/neuonc/nou138. Epub 2014 Jul 15. PubMed
PMID: 25028501; PubMed Central PMCID: PMC4232086.

  * Paulsson AK, Holmes JA, Peiffer AM, Miller LD, Liu W, Xu J, Hinson WH, Lesser 
GJ, Laxton AW, Tatter SB, Debinski W, Chan MD. Comparison of clinical outcomes
and genomic characteristics of single focus and multifocal glioblastoma. J
Neurooncol. 2014 Sep;119(2):429-35. doi: 10.1007/s11060-014-1515-1. Epub 2014 Jul
3. PubMed PMID: 24990827; PubMed Central PMCID: PMC4146694.

  * Molenaar RJ, Verbaan D, Lamba S, Zanon C, Jeuken JW, Boots-Sprenger SH,
Wesseling P, Hulsebos TJ, Troost D, van Tilborg AA, Leenstra S, Vandertop WP,
Bardelli A, van Noorden CJ, Bleeker FE. The combination of IDH1 mutations and
MGMT methylation status predicts survival in glioblastoma better than either IDH1
or MGMT alone. Neuro Oncol. 2014 Sep;16(9):1263-73. doi: 10.1093/neuonc/nou005.
Epub 2014 Feb 6. PubMed PMID: 24510240; PubMed Central PMCID: PMC4136888.

  * Kanemoto M, Shirahata M, Nakauma A, Nakanishi K, Taniguchi K, Kukita Y,
Arakawa Y, Miyamoto S, Kato K. Prognostic prediction of glioblastoma by
quantitative assessment of the methylation status of the entire MGMT promoter
region. BMC Cancer. 2014 Aug 30;14:641. doi: 10.1186/1471-2407-14-641. PubMed
PMID: 25175833; PubMed Central PMCID: PMC4161852.

  * Collins VP, Ichimura K, Di Y, Pearson D, Chan R, Thompson LC, Gabe R, Brada
M, Stenning SP; BR12 Collaborators. Prognostic and predictive markers in
recurrent high grade glioma; results from the BR12 randomised trial. Acta
Neuropathol Commun. 2014 Jun 20;2:68. doi: 10.1186/2051-5960-2-68. PubMed PMID:
24952577; PubMed Central PMCID: PMC4229733.

  * Xu M, Nekhayeva I, Cross CE, Rondelli CM, Wickliffe JK, Abdel-Rahman SZ.
Influence of promoter/enhancer region haplotypes on MGMT transcriptional
regulation: a potential biomarker for human sensitivity to alkylating agents.
Carcinogenesis. 2014 Mar;35(3):564-71. doi: 10.1093/carcin/bgt355. Epub 2013 Oct 
25. PubMed PMID: 24163400; PubMed Central PMCID: PMC3941746.

  * Quillien V, Lavenu A, Sanson M, Legrain M, Dubus P, Karayan-Tapon L, Mosser
J, Ichimura K, Figarella-Branger D. Outcome-based determination of optimal
pyrosequencing assay for MGMT methylation detection in glioblastoma patients. J
Neurooncol. 2014 Feb;116(3):487-96. doi: 10.1007/s11060-013-1332-y. Epub 2014 Jan
14. PubMed PMID: 24420923; PubMed Central PMCID: PMC3905192.

  * Hiddingh L, Tannous BA, Teng J, Tops B, Jeuken J, Hulleman E, Boots-Sprenger 
SH, Vandertop WP, Noske DP, Kaspers GJ, Wesseling P, Wurdinger T. EFEMP1 induces 
γ-secretase/Notch-mediated temozolomide resistance in glioblastoma. Oncotarget.
2014 Jan 30;5(2):363-74. PubMed PMID: 24495907; PubMed Central PMCID: PMC3964213.

  * Shen D, Liu T, Lin Q, Lu X, Wang Q, Lin F, Mao W. MGMT promoter methylation
correlates with an overall survival benefit in Chinese high-grade glioblastoma
patients treated with radiotherapy and alkylating agent-based chemotherapy: a
single-institution study. PLoS One. 2014 Sep 11;9(9):e107558. doi:
10.1371/journal.pone.0107558. eCollection 2014. PubMed PMID: 25211033; PubMed
Central PMCID: PMC4161443.

  * Lai RK, Chen Y, Guan X, Nousome D, Sharma C, Canoll P, Bruce J, Sloan AE,
Cortes E, Vonsattel JP, Su T, Delgado-Cruzata L, Gurvich I, Santella RM, Ostrom
Q, Lee A, Gregersen P, Barnholtz-Sloan J. Genome-wide methylation analyses in
glioblastoma multiforme. PLoS One. 2014 Feb 21;9(2):e89376. doi:
10.1371/journal.pone.0089376. eCollection 2014. PubMed PMID: 24586730; PubMed
Central PMCID: PMC3931727.

  * Oberstadt MC, Bien-Möller S, Weitmann K, Herzog S, Hentschel K, Rimmbach C,
Vogelgesang S, Balz E, Fink M, Michael H, Zeden JP, Bruckmüller H, Werk AN,
Cascorbi I, Hoffmann W, Rosskopf D, Schroeder HW, Kroemer HK. Epigenetic
modulation of the drug resistance genes MGMT, ABCB1 and ABCG2 in glioblastoma
multiforme. BMC Cancer. 2013 Dec 31;13:617. doi: 10.1186/1471-2407-13-617. PubMed
PMID: 24380367; PubMed Central PMCID: PMC3890604.

  * Brennan CW, Verhaak RG, McKenna A, Campos B, Noushmehr H, Salama SR, Zheng S,
Chakravarty D, Sanborn JZ, Berman SH, Beroukhim R, Bernard B, Wu CJ, Genovese G, 
Shmulevich I, Barnholtz-Sloan J, Zou L, Vegesna R, Shukla SA, Ciriello G, Yung
WK, Zhang W, Sougnez C, Mikkelsen T, Aldape K, Bigner DD, Van Meir EG, Prados M, 
Sloan A, Black KL, Eschbacher J, Finocchiaro G, Friedman W, Andrews DW, Guha A,
Iacocca M, O'Neill BP, Foltz G, Myers J, Weisenberger DJ, Penny R, Kucherlapati
R, Perou CM, Hayes DN, Gibbs R, Marra M, Mills GB, Lander E, Spellman P, Wilson
R, Sander C, Weinstein J, Meyerson M, Gabriel S, Laird PW, Haussler D, Getz G,
Chin L; TCGA Research Network. The somatic genomic landscape of glioblastoma.
Cell. 2013 Oct 10;155(2):462-77. doi: 10.1016/j.cell.2013.09.034. Erratum in:
Cell. 2014 Apr 24;157(3):753. PubMed PMID: 24120142; PubMed Central PMCID:
PMC3910500.

  * Thon N, Kreth S, Kreth FW. Personalized treatment strategies in glioblastoma:
MGMT promoter methylation status. Onco Targets Ther. 2013 Sep 27;6:1363-72. doi: 
10.2147/OTT.S50208. Review. PubMed PMID: 24109190; PubMed Central PMCID:
PMC3792931.

  * Wiestler B, Claus R, Hartlieb SA, Schliesser MG, Weiss EK, Hielscher T,
Platten M, Dittmann LM, Meisner C, Felsberg J, Happold C, Simon M, Nikkhah G,
Papsdorf K, Steinbach JP, Sabel M, Grimm C, Weichenhan D, Tews B, Reifenberger G,
Capper D, Müller W, Plass C, Weller M, Wick W; Neuro-oncology Working Group (NOA)
of the German Cancer Society. Malignant astrocytomas of elderly patients lack
favorable molecular markers: an analysis of the NOA-08 study collective. Neuro
Oncol. 2013 Aug;15(8):1017-26. doi: 10.1093/neuonc/not043. Epub 2013 Apr 17.
PubMed PMID: 23595628; PubMed Central PMCID: PMC3714152.

  * Qiu S, Lin S, Hu D, Feng Y, Tan Y, Peng Y. Interactions of
miR-323/miR-326/miR-329 and miR-130a/miR-155/miR-210 as prognostic indicators for
clinical outcome of glioblastoma patients. J Transl Med. 2013 Jan 9;11:10. doi:
10.1186/1479-5876-11-10. PubMed PMID: 23302469; PubMed Central PMCID: PMC3551827.

  * Weller M, Stupp R, Hegi ME, van den Bent M, Tonn JC, Sanson M, Wick W,
Reifenberger G. Personalized care in neuro-oncology coming of age: why we need
MGMT and 1p/19q testing for malignant glioma patients in clinical practice. Neuro
Oncol. 2012 Sep;14 Suppl 4:iv100-8. doi: 10.1093/neuonc/nos206. Review. PubMed
PMID: 23095825; PubMed Central PMCID: PMC3480248.

  * McDonald KL, Aw G, Kleihues P. Role of Biomarkers in the Clinical Management 
of Glioblastomas: What are the Barriers and How Can We Overcome Them? Front
Neurol. 2013 Jan 18;3:188. doi: 10.3389/fneur.2012.00188. eCollection 2012.
PubMed PMID: 23346075; PubMed Central PMCID: PMC3548232.

  * Zinn PO, Sathyan P, Mahajan B, Bruyere J, Hegi M, Majumder S, Colen RR. A
novel volume-age-KPS (VAK) glioblastoma classification identifies a prognostic
cognate microRNA-gene signature. PLoS One. 2012;7(8):e41522. doi:
10.1371/journal.pone.0041522. Epub 2012 Aug 3. PubMed PMID: 22870228; PubMed
Central PMCID: PMC3411674.


