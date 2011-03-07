\name{artificial.data
}
\alias{artificial.data
}
\title{Create an artificial data set
}
\description{artificial.data creates an artifical data set of multispecies presence absence data based on a binomial response
}
\usage{artificial.data(formula,data,theta,S)
}
\arguments{
\item{ formula}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted}
\item{ data}{a data frame containing the variables in the formula}
\item{ theta}{coefficents for each species archetype. Matrix of G x number of parameters. Each row is a different species archetype.}
\item{ S}{number of species to create}
}
\details{each species is a random realisation of the archetype species
}
\value{a list with elements pa (matrix of presence absence data) and group (membership of each species to each group)
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
## for bernoulli
theta <- matrix(c(-0.9,-0.6,0.5,1,-0.9,1,0.9,-0.9),4,2,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x=runif(100,0,2.5))
dat1 <- artificial.data(y~1+x,dat,theta,20)
## for negbin
theta <- matrix(c(-0.9,-0.6,0.5,1,-0.9,1,0.9,-0.9),4,2,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x=runif(100,0,2.5))
dat1 <- artificial.data(y~1+x,dat,theta,20)
}
\keyword{misc}
