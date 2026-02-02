\name{permalgorithm}
\alias{permalgorithm}
\title{Generate Event Times Conditional On Time-Dependent Covariates}
\description{This version of the permutational algorithm generates a dataset in which event and censoring times are conditional on an user-specified list of covariates, some or all of which are time-dependent.  Event times and censoring times also follow user-specified distributions.}
\usage{permalgorithm(numSubjects, maxTime, Xmat, XmatNames = NULL,
eventRandom = NULL, censorRandom = NULL, betas, groupByD = FALSE)}
\arguments{
  \item{numSubjects}{is the number of subjects generated.}
  \item{maxTime}{is a non-zero integer represening the maximum length
    of follow-up.}
  \item{Xmat}{is the matrix of covariates values in a counting process
  format where every line represent one and only one time interval,
  during which all covariate values for a given subject remains
  constant. Consequently, \code{Xmat} should have
  \code{numSubjects}*\code{maxTime} rows. Each column of \code{Xmat}
  corresponds to a different covariates on which the event is conditionned. For fixed-in-time covariates, the same value should be
  replicated in each of \code{maxTime} row for a given subject.}
  \item{XmatNames}{a an optional vector of character strings representing
  the names of each of the covariates in \code{Xmat}.}
  \item{eventRandom}{represents individual event times. \code{eventRandom}
  can be a vector of nonegative integer values or a random generating
  function with argument n. In both cases, its values must be smaller or
  equal to \code{maxTime}. If left unspecified, then the algorithm
  generates event times based on an uniform distribtuion [1, maxTime].}
  \item{censorRandom}{represents individual censoring
  times. \code{censorRandom} can be a vector of nonegative integer
  values or a random generating function with argument n. In both cases,
  its values must be smaller or equal to \code{maxTime}.  The default is
  Uniform[1,maxTime].}
  \item{betas}{is a vector of regression coefficients (log hazard) that represent the
  magnitude of the relationship between each of the covariates and the
  risk of an event. The length of \code{betas} should correspond to the
  number of columns in \code{Xmat}.}
  \item{groupByD}{groupByD is an option that, when enabled, increases the
  computational efficiency of the algorithm by replacing the individual
  assignment of event times and censoring times by grouped
  assignements. The side effect of this option is that it generates
  datasets that are, on average, slightly less consistent with the model
  described by \code{betas} that those generated  with the groupByD
  option set to FALSE. Still, groupByD=TRUE may be  usefull to generate
  large datasets where \code{maxTime} is much smaller than
  \code{numSubjects} so that many ties are expected. Default is FALSE.}
}
\details{
The gist of the algorithm is to perform a one-to-one matching of
\code{n} observed times with independently generated vectors of
covariates values.  The matching is performed based on a permutation
probability law derived from the partial likelihood of Cox's 
Proportional Hazards (PH) model.  

The number of events obtained in the data.frame returned by the function
depends on both the distribution of event \code{enventRandom} and
censoring times \code{censorRandom}.  In the simplest case where the
distirbution of \code{eventRandom} is Uniform over follow-up U[1,m], and
the censoring is random, the number of observed events in the data.frame
returnd by the algorithm is determined by the upper bound of the Uniform
distribution of \code{censorRandom}.  For example, setting the
distribution of \code{censorRandom} to U[1,m] will lead to approximately
half of the subjects to experience an event during follow-up, while
setting the distribution of \code{censorRandom} to U[1,3/2] will lead to
approximately two thirds of the observed times to be events.

Subjects without an event before or on \code{maxTime} and who are not
censored before \code{maxTime} are censored on \code{maxTime}
(administrative censoring).

*** Warning *** Currently the algorithm only takes Xmat in matrix
format. Consequently, factor variables are not allowed. Instead, users need
to code them with binary indicators. 

}
\value{
  A data.frame object with columns corresponding to
  \item{Id}{Identifies the rows of the data.frame that corresponds to
    each of the \code{n} individuals.}
  \item{Event}{Indicator of event.  \code{Event} = 1 when event occurs
    and 0 otherwise.}
  \item{Fup}{Individual follow-up time.}
  \item{Start}{For counting process formulation.  Represents the start
    of each time interval.}
  \item{Stop}{For counting process formulation.  Represents the end of
    each time interval.}
  \item{Xmat}{The values of the covariates specified in  Xmat.}
}

\references{This algorithm is an extension of the permutational algorithm 
first introduced by Abrahamowicz, MacKenzie and Esdaile, and described in
details by MacKenzie and Abrahamowicz.  The current version of the permutational
algorithm is a flexible tool to generate event and censoring times that follow
user-specified distributions and that are conditional on user-specified
covariates. This is especially useful whenever at least one of the
covariate is time-dependent so that conventional inversion methods are
difficult to implement.

The algorithm has been validated through simulations in Sylvestre and
Abrahamowicz.  Please reference the manuscript by Sylvestre and Abrahamowicz,
cited below, if this program is used in any published material.

Sylvestre M.-P., Abrahamowicz M. (2008) Comparison of algorithms to generate
event times conditional on time-dependent covariates. \emph{Statistics in 
  Medicine} \bold{27(14)}:2618--34.

Abrahamowicz M., MacKenzie T., Esdaile J.M. (1996) Time-dependent hazard ratio: 
modelling and hypothesis testing with application in lupus nephritis.
\emph{JASA} \bold{91}:1432--9.

MacKenzie T., Abrahamowicz M. (2002) Marginal and hazard ratio specific random
 data generation: Applications to semi-parametric bootstrapping.
\emph{Statistics and Computing} \bold{12(3)}:245--252.

 }


\author{Marie-Pierre Sylvestre, Thad Evans, Todd MacKenzie, Michal Abrahamowicz}
\examples{

# Example 1 - Generating adverse event conditional on use
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


# Example 2 - Generating Myocardial Infarction (MI) conditional on
# biennial measures of systolic blood pressure (like in the
# Framingham data).

m = 16 # exams
n <- 10000 # individuals

# Very crude way to generate the data,  meant as an example only!
sysBP <- rnorm(n*m,  120, 15) 

# by not submitting event and censor time, one let the algorithm
# generate them from uniform distributions over the follow-up time.

data2 <- permalgorithm(n, m, sysBP, XmatNames="sysBP", betas=log(1.01),
groupByD=FALSE )

}

\keyword{survival}

