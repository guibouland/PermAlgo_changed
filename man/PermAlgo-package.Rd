\name{PermAlgo-package}
\alias{PermAlgo-package}
\alias{PermAlgo}
\docType{package}
\title{Generate Event Times Conditional On Time-Dependent Covariates
}
\description{
This version of the permutational algorithm generates a dataset in which
event and censoring times are conditional on an user-specified list of
covariates, some or all of which are time-dependent.  Event times and censoring
times also follow user-specified distributions. 
}
\details{
\tabular{ll}{
Package: \tab PermAlgo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-08-24\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package contains one function avialable to the user, \code{permalgorithm}. The gist of the algorithm is to perform a one-to-one matching of \code{n} observed times with \code{n} independently generated vectors of covariates values.  The matching is performed based on a permutation probability law derived from the partial likelihood of Cox's Proportional Hazards (PH) model.  
}
\author{Marie-Pierre Sylvestre, Thad Evans, Todd MacKenzie, Michal Abrahamowicz

Maintainer: Marie-Pierre Sylvestre  <marie-pierre.sylvestre.chum@ssss.gouv.qc.ca>
}
\references{
This algorithm is an extension of the permutational algorithm 
first introduced by Abrahamowicz, MacKenzie and Esdaile, and described in
details by MacKenzie and Abrahamowicz.  The current version of the permutational
algorithm is a flexible tool to generate event and censoring times that follow
user-specified distributions and that are conditional on user-specified
covariates.  It has been validated through simulations in Sylvestre and
Abrahamowicz.  Please reference the manuscript by Sylvestre and Abrahamowicz
cited below if the results of this program are used in any published
material.

Sylvestre M.-P., Abrahamowicz M. (2008) Comparison of algorithms to generate
event times conditional on time-dependent covariates. \emph{Statistics in
  Medicine} \bold{27(14)}:2618--34.

Abrahamowicz M., MacKenzie T., Esdaile J.M. (1996) Time-dependent hazard ratio:
modelling and hypothesis testing with application in lupus nephritis. \emph{JASA} 
\bold{91}:1432--9.

MacKenzie T., Abrahamowicz M. (2002) Marginal and hazard ratio specific random
data
generation: Applications to semi-parametric bootstrapping. \emph{Statistics
and Computing} \bold{12(3)}:245--252.
}

\examples{

# Example - Generating adverse event conditional on use
# of prescription drugs

# Prepare the matrice of covariate (Xmat)
# Here we simulate daily exposures to 2 prescription drugs over a
# year. Drug prescriptions can start any day of follow-up, and their
# duration is a multiple of 7 days. There can be multiple prescriptions
# for each individuals over the year and interuptions of drug use in
# between.

# Additionaly, there is a time-independant binary covarite (sex).

n=500 # subjects
m=365 # days

# Generate the matrix of three covariate, in a 'long' format.
Xmat=matrix(ncol=3, nrow=n*m)

# time-independant binary covariate
Xmat[,1] <- rep(rbinom(n, 1, 0.3), each=m)

# Function to generate an individual time-dependent exposure history
# e.g. generate prescriptions of different durations and doses.
TDhist <- function(m){
  start <- round(runif(1,1,m),0) # individual start date
  duration <-  7 + 7*rpois(1,3) # in weeks
  dose <-  round(runif(1,0,10),1) 
  vec <- c(rep(0, start-1), rep(dose, duration))
  while (length(vec)<=m){
    intermission <- 21 + 7*rpois(1,3) # in weeks
    duration <-  7 + 7*rpois(1,3) # in weeks
    dose <-  round(runif(1,0,10),1)
    vec <- append(vec, c(rep(0, intermission), rep(dose, duration)))}
  return(vec[1:m])}

# create TD var
Xmat[,2] <- do.call("c", lapply(1:n, function(i) TDhist(m)))
Xmat[,3] <- do.call("c", lapply(1:n, function(i) TDhist(m)))

# genereate vectors of event and censoring times prior to calling the
# function for the algorithm

eventRandom <- round(rexp(n, 0.012)+1,0)
censorRandom <- round(runif(n, 1,870),0)

# Generate the survival data conditional on the three covariates

data <- permalgorithm(n, m, Xmat, XmatNames=c("sex", "Drug1", "Drug2"),
eventRandom = eventRandom, censorRandom=censorRandom, betas=c(log(2),
log(1.04), log(0.99)), groupByD=FALSE )

# could use survival library and check whether the data was generated
# properly using coxph(Surv(Start, Stop, Event) ~ sex + Drug1 + Drug2,
# data)

}
\keyword{ survival }
