\name{print.archetype}
\alias{print.archetype}
\title{Print archetypes from fitted model}
\description{print.archetype
}
\section{Method}{print(x,..)
}
\arguments{
\item{ x }{a model object returned from SpeciesMix}
\item{ ... }{further arguments passed to or from other methods.}
}
\details{Prints the fitted model returned from SpeciesMix, giving the coefficents, standard errors and posterior probabilities
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
print(fm4)
}
\keyword{misc}
