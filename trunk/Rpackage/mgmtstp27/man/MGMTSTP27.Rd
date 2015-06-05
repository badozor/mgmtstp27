\name{MGMTSTP27}
\alias{MGMTSTP27}
\docType{data}
\title{
Model to predict MGMTSTP27
}
\description{

}
\usage{data("MGMTSTP27")}
\format{
  an object 'glm' with additional numerical vector called 'perf' containing performance information.
  \describe{
    \item{\code{perf}}{a numeric vector containing information}
  }
}
\details{
The function MGMTpredict provides prediction of DNA methylation status of
MGMT promoter as described in Bady et al (2012). The model is an object glm with an additional numerical vector called 'perf' containing performance information.}
\source{
The model is presented in details in Bady et al. (2012).
}
\references{
Bady, P., D. Sciuscio, A.-C. Diserens, J. Bloch, M. J. van den Bent, C. Marosi, P.-Y. Dietrich, M. Weller, L. Mariani, F. L. Heppner, D. R. McDonald, D. Lacombe, R. Stupp, M. Delorenzi, and M. E. Hegi. 2012. MGMT methylation analysis of glioblastoma on the Infinium methylation BeadChip identifies two distinct CpG regions associated with gene silencing and outcome, yielding a prediction model for comparisons across datasets, tumor grades, and CIMP-status. Acta Neuropathologica 124:547-560.
}
\examples{
data(MGMTSTP27)
MGMTSTP27
# maximization of good classification
MGMTSTP27$perf1
# balance amon specificity and sensitivity
MGMTSTP27$perf2
}
\keyword{datasets}
