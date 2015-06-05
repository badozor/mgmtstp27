\name{mgmt}
\alias{MGMTpredict}
\alias{MGMTsim}
\alias{MGMTqc}
\alias{MGMTqc.pop}
\alias{MGMTqc.single}
\title{
set of tools related to prediction of the DNA methylation of MGMT promoter.
}
\description{
set of tools related to prediction of the DNA methylation of MGMT promoter.
}
\usage{
MGMTpredict(x,level = 0.05,dispersion=FALSE,transpose=FALSE,ic.distrib="normal",cutoff=1,...)
MGMTsim(n=1000,proba=NULL,newdata=NULL,...)
MGMTqc.pop(object,sim=FALSE,n=1000,which.plot=1:6,mfrow=c(3,3),...)
MGMTqc.single(object,nsample=NULL,sim=FALSE,n=1000,which.plot=1:4,mfrow=c(2,3),...)
}
\arguments{
  \item{x}{
a data.frame containing contining the M-values for the probes 'cg12981137'  and 'cg12434587' (in row and columns).
}
  \item{object}{
an object of the class 'mgmt'.
}
  \item{level}{
a numeric value corresponding to level used to compute confidence intervals (level=0.05, by default).
}
  \item{dispersion}{
a logical value. If dispersion=TRUE, dispersion correction was used to compute Confidence intervals.
}
  \item{transpose}{
a logical value (transpose=TRUE, when the probes are organized by rows).
}
  \item{ic.distrib}{
a character value indicating the probability distribution used to compute confidence intervals 
("student" or "normal", by default ic.distrib="normal").
}
  \item{cutt-off}{
a numeric value identifying the cut-off used to calculate the DNA methylation state of MGMT promoter 
(1= better classification, 2= balance among sensitivity and specificty). By default, cutoff is equal to 1.
}
  \item{n}{
number of simulated samples
}
%  \item{proba}{
%prevalance of methylated MGMT promoter.
%}
  \item{newdata}{
an object of the class 'mgmt' containing containing the M-values for the probes 'cg12981137'  and 'cg12434587' (in columns).
}
  \item{sim}{
a logical value.If sim=TRUE, the reference is obtained by simulation (Multivariate distribution. If sim=FALSE (by default), 
the reference is the training data.
}
  \item{nsample}{
a numerical value identifying the sample of interest (only used in the function MGMTqc.single, bydefault nsample=1).
}
  \item{which.plot}{
a vector telling which plots to produce.
}
  \item{mfrow}{
parameter for the array of figures to be drawn .
}
  \item{...}{
further arguments passed to or from other methods
}

}
\details{
information about simulation and QC grphical ouput
}
\value{
The function 'MGMTsim' return a data.frame containing two simulated M-values corresponding to the probes 'cg12981137'  
and 'cg12434587'. The function 'MGMTpredict' return data.frame contianing the following information:
  \item{sample}{sample(row) names}
  \item{cg12434587}{M-value for the probe 'cg12434587'}
  \item{cg12981137}{M-value for the probe 'cg12981137'}
  \item{pred}{probability that MGMT promoter is methylated}
  \item{lower}{lower limit of the confidence intervals for the probability}
  \item{upper}{upper limit of the confidence intervals for the probability}
  \item{state}{DNA methylation state of MGMT promoter using the cut-off provided in Bady et al. (2012)}  
}
\references{
Bady, P., D. Sciuscio, A.-C. Diserens, J. Bloch, M. J. van den Bent, C. Marosi, P.-Y. Dietrich, M. Weller, L. Mariani, 
F. L. Heppner, D. R. McDonald, D. Lacombe, R. Stupp, M. Delorenzi, and M. E. Hegi. 2012. 
MGMT methylation analysis of glioblastoma on the Infinium methylation BeadChip identifies two distinct 
CpG regions associated with gene silencing and outcome, yielding a prediction model for comparisons across datasets, 
tumor grades, and CIMP-status. Acta Neuropathologica 124:547-560.
}
\author{
P.BADY
}
\examples{
data(MGMTSTP27)
training1 <- MGMTSTP27$data
pred1 <- MGMTpredict(training1)
sim1 <- MGMTsim(n=100,newdata=pred1) 
qqplot(pred1[,"cg12434587"],sim1[,"cg12434587"])
MGMTqc(pred1)
}
\keyword{mgmt}
