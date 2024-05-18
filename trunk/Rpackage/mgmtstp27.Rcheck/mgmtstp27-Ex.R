pkgname <- "mgmtstp27"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('mgmtstp27')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MGMTSTP27")
### * MGMTSTP27

flush(stderr()); flush(stdout())

### Name: MGMTSTP27
### Title: Model to predict MGMTSTP27
### Aliases: MGMTSTP27
### Keywords: datasets

### ** Examples

data(MGMTSTP27)
MGMTSTP27
# maximization of good classification
MGMTSTP27$perf1
# balance amon specificity and sensitivity
MGMTSTP27$perf2



cleanEx()
nameEx("gbm")
### * gbm

flush(stderr()); flush(stdout())

### Name: gbm
### Title: DNA methylation of MGMT promoter region from Infinium HM-450k
###   and HM-27K platforms
### Aliases: gbm NCHgbm450 TCGAgbm27
### Keywords: datasets

### ** Examples

# table S3 (bady et al 2012)
data(NCHgbm450)
head(NCHgbm450)
# table S4 (bady et al 2012)
data(TCGAgbm27)
head(TCGAgbm27)



cleanEx()
nameEx("mgmt")
### * mgmt

flush(stderr()); flush(stdout())

### Name: mgmt
### Title: set of tools related to prediction of the DNA methylation of
###   MGMT promoter.
### Aliases: MGMTpredict MGMTsim MGMTqc MGMTqc.pop MGMTqc.single
### Keywords: mgmt

### ** Examples

data(MGMTSTP27)
training1 <- MGMTSTP27$data
pred1 <- MGMTpredict(training1)
sim1 <- MGMTsim(n=100,newdata=pred1) 
qqplot(pred1[,"cg12434587"],sim1[,"cg12434587"])
MGMTqc.pop(pred1)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
