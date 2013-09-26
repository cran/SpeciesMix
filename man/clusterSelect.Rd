\name{clusterSelect}
\alias{clusterSelect}
\title{Fit SpeciesMix for a range of vaules for G}
\description{clusterSelect fits models with varying values of G to determine the appropriate number of archetype species.
}
\usage{
#clusterSelect(sp.form,sp.data,covar.data,G=1:10,\\n
#em.prefit=TRUE, em.steps=4 ,em.refit=3,\\n
#est.var=FALSE,trace=TRUE)
}
\arguments{
\item{ sp.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted}
\item{ sp.data}{a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned.}
\item{ covar.data}{a data frame containng the covariate data for each site. Names of columns must match that given in}
\item{ G}{Vector containing the range of archetype species to fit.}
\item{ em.prefit}{obtain initial parameter estimates from EM}
\item{ em.steps}{number of EM steps to do if using em.prefit}
\item{ em.refit}{refits model so that the global maxima can be found using EM.}
\item{ est.var}{calculate the variance covariace matrix for each group}
\item{ trace}{the trace of the EM steps}
}
\details{fits multiple fitMix models across the range of values for G. Most of the arguments are passed directly to fitMix
}
\value{
\item{ aic}{vector containing the aic value for each value of G}
\item{ bic}{bic}
\item{ fm}{a list containing all output from each vaule of G.}
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
G <-4
S <- 20
theta <- matrix(c(-0.9,-0.6,0.5,1,-0.9,1,0.9,-0.9),4,2,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x=runif(100,0,2.5),z=rnorm(100,10,2))
dat1 <- artificial.data(y~1+x,dat,theta,S)
dat <- dat[,2:3]
clusters <- clusterSelect(obs~1+x,dat1$pa,dat,G=2:5,em.refit=2)
}
\keyword{misc}
