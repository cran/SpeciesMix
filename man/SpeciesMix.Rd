\name{SpeciesMix
}
\alias{SpeciesMix
}
\title{Fit finite mixture models to species distributions
}
\description{SpeciesMix is used to fit mixtures of glms to multivariate species data. The function uses BFGS to optimise the mixture likelihood.There is the option to use EM get appropriate starting parameters.
}
\usage{SpeciesMix(sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=4, em.refit = 1 , est.var = FALSE,residuals=FALSE,trace=TRUE)
}
\arguments{
\item{ sp.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted. The response variable is arbitary, all the responses are in sp.data.}
\item{ sp.data}{a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned.}
\item{ covar.data}{a data frame containng the covariate data for each site. Names of columns must match that given in sp.form}
\item{ G}{number of archetype species to fit}
\item{ pars}{vector of parameters order is additive logisic transform of pi, glm coefficents by column (ie all intercepts first, etc)}
\item{ em.prefit }{find initial parameter values using EM}
\item{ em.steps }{number of EM steps to preform}
\item{ em.refit }{number of times to refit using EM}
\item{ est.var }{estimate the variance covariance matrix (TRUE/FALSE)}
\item{ residuals }{claculate species wise residuals}
\item{ trace}{the trace}
}
\details{SpeciesMix is a wrapper for fitmix.cpp that allows for easier data input. The data frames are merged into the appropriate format for the use in fitmix.cpp. Minima is found using vmmin (BFGS) and the gradients are calculated using CPPAD (auto differentiation)
}
\value{
\item{ logl}{loglikelihood}
\item{ aic}{AIC}
\item{ bic}{BIC}
\item{ tau }{posterior probabilities giving group membership. These values indicate which archetype each species belongs to.}
\item{ pi }{prior probabilities of archetype membership.}
\item{ coef}{coefficents for each species archetype}
\item{ covar}{variance co-variance matrix at the maxima. Calculated if est.var is TRUE}
\item{ hessian}{Hessian matrix at maximum likelihood. Calculated if est.var is TRUE}
\item{ gradient}{derivitive of parameters at the maxima. Calculated if est.var is TRUE}
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
G <-4
S <- 50
theta <- matrix(c(-9,35,-32,0,0.7,0,-16,23,-8.2,-3,-0.6,0.8),4,3,byrow=TRUE)
dat <- data.frame(y=rep(1,200),x=runif(200,0,2.5),z=rnorm(200,10,2))
dat <- data.frame(dat,x.sq=dat$x^2)
dat1 <- artificial.data(y~1+x+x.sq,dat,theta,S)
fm4 <- SpeciesMix(obs~1+x+x.sq,dat1$pa,dat,G=4,em.prefit=TRUE,em.refit=1,est.var=TRUE)
}
\keyword{misc}
