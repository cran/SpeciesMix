\name{logLike.pars}
\alias{logLike.pars}
\title{logLike.pars}
\description{logLike.pars return loglike for give parameters. pi is on normal scale with G groups}
\usage{
logLike.pars(pi,coef,sp.form,sp.data,covar.data)
}
\arguments{
\item{ pi}{vector of pi returned from SpeciesMix}
\item{ coef}{matrix of coefficents returned from Species Mix}
\item{ sp.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted}
\item{ sp.data}{a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned.}
\item{ covar.data}{a data frame containng the covariate data for each site. Names of columns must match that given in}
}
\details{To Come}
\value{\item{ logl}{loglikelihood}}
\section{Authors}{Piers Dunstan and Scott Foster}
\examples{
G <-4
S <- 50
theta <- matrix(c(-9,35,-32,0,0.7,0,-16,23,-8.2,-3,-0.6,0.8),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,200),x=runif(200,0,2.5),z=rnorm(200,10,2))
dat <- data.frame(dat,x.sq=dat$x^2)
dat1 <- artificial.data(y~1+x+x.sq,dat,theta,S)
fm4 <- SpeciesMix(obs~1+x+x.sq,dat1$pa,dat,G=4,em.prefit=TRUE,em.refit=1,est.var=FALSE)
logLike.pars(fm4$pi,fm4$coef,obs~1+x+x.sq,dat1$pa,dat)
}
\keyword{misc}
