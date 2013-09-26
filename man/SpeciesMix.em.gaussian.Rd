\name{SpeciesMix.em.gaussian}
\alias{SpeciesMix.em.gaussian}
\title{Implementation of SpeciesMix using EM algorithim}
\description{
SpeciesMix.em is used to fit mixtures of glms to multivariate species data. The function is a wrapper function for fitMix but allows for multiple EM starts to find the global maxima.

Do not use this function under normal conditions. Documentation is included to allow direct interface with the EM steps.
}
\usage{
#SpeciesMix.em.gaussian(sp.form,sp.data,covar.data,G=2,\\n
#em.refit=1,ite.max=500, est.var = FALSE,\\n
#residuals=FALSE,trace=TRUE)
}
\arguments{
\item{ sp.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted}
\item{ sp.data}{a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned.}
\item{ covar.data}{a data frame containng the covariate data for each site. Names of columns must match that given in}
\item{ G}{number of archetype species to fit}
\item{ em.refit}{refits model so that the global maxima can be found using EM.}
\item{ ite.max}{maximum number of EM steps to preform}
\item{ est.var}{calculate variance covariance matrix}
\item{ residuals}{calculate residuals}
\item{ trace}{the trace of the EM steps}
}
\details{SpeciesMix is a wrapper for fitMix that allows for easier data input. The data frames are merged into the appropriate format for the use in fitMix. Multiple models can be fitted so that
}
\value{
\item{ logl}{loglikelihood}
\item{ aic}{AIC}
\item{ bic}{BIC}
\item{ ICL}{Integrated Classification Likelihood Criterion}
\item{ tau}{posterior probabilities giving group membership}
\item{ coef}{coefficents for each species archetype}
\item{ fmM}{list of glm objects for archetypes 1 to G}
\item{ model.tau}{full matrix of taus used by glm for model fitting. May be useful for residuals}
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
G <-4
S <- 20
theta <- matrix(c(-0.9,-0.6,0.5,1,-0.9,1,0.9,-0.9),4,2,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x=runif(100,0,2.5),z=rnorm(100,10,2))
dat1 <- artificial.data(y~1+x,dat,theta,S)
}
\keyword{misc}
