\name{predict.archetype}
\alias{predict.archetype}
\title{Predict presence of archetypes from fitted model}
\description{predict.archetype
}
\section{Method}{predict(object, new.obs,\dots{})
}
\arguments{
\item{ object }{a model object returned from SpeciesMix}
\item{ new.obs}{covariate data at new locations for predictions. Must be a data frame with same row names as in the model formula}
\item{ ... }{further arguments passed to or from other methods.}
}
\details{predicts probability of occurance at new location with standard errors calculated using delta method
}
\value{
\item{ fit}{predicted probability of presence}
\item{ se.fit}{standard error of fit}
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
new.obs <- data.frame(x=seq(0,2.5,0.01),x.sq=seq(0,2.5,0.01)^2)
predict(fm4,new.obs)
}
\keyword{misc}
