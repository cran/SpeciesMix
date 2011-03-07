# This is package SpeciesMix 

"additive.logistic" <-
function (x,inv=FALSE) 
{
  if(inv){
    x <- log(x/x[length(x)])
    return(x)
  }

  x.t <- exp(x)
  x.t <- x.t/(1+sum(x.t))
  x.t[length(x.t)+1] <- 1-sum(x.t)
  return(x.t)
}


"apply.glm" <-
function (i,form,datsp,tau,n) 
{
  dat.tau <- rep(tau[,i],each=n)
  x <- model.matrix(as.formula(form),data=datsp)
  y <- datsp$obs
  ##dat.tau <- get("dat.tau")
  ##datsp <- get("datsp")
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  f.mix <- glm.fit(x,y,dat.tau,family=binomial())
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
  ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
  list(coef=f.mix$coef)
}


"artificial.data" <-
function (formula,data,theta,S) 
{
  X <- model.matrix(formula,data)
  out <- matrix(0,dim(X)[1],S)
  k <- dim(theta)[1] ## rows of theta = number of groups
  group <- rep(0,S)
  for(s in 1:S){
    g <- ceiling(runif(1)*k)
    lgtp <- X%*%theta[g,]
    p <- exp(lgtp)/(1+exp(lgtp))
    out[,s] <- rbinom(dim(X)[1],1,p)
    group[s] <- g
  }
  pi <- tapply(group,group,length)/S
  list(pa=out,group=group,pi=pi) 
}




"c.f.n" <-
function(x){as.numeric(as.character(x))}


"clusterSelect" <-
function (sp.form,sp.data,covar.data,G=1:10,em.prefit=TRUE, em.steps=4 ,em.refit=3,est.var=FALSE,trace=TRUE) 
{
  my.fun <- function(g,form,sp.data,covar.data){
    cat("Fitting group",g,"\n")
    try(SpeciesMix(form,sp.data,covar.data,g,em.prefit=em.prefit,em.steps=em.steps,em.refit=em.refit,est.var=est.var,trace=trace))
  }
  

#  if(mc){ out <- mclapply(G,my.fun,form,dat,mc.preschedule = FALSE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = set.cores)} else
  { out <- lapply(G,my.fun,sp.form,sp.data,covar.data)}

  aic <- rep(0,length(G))
  bic <- rep(0,length(G))
  fm <- list()
  for(i in 1:length(G))
    if(!is.atomic(out[[i]])){
      aic[i] <- out[[i]]$aic
      bic[i] <- out[[i]]$bic
      fm[[i]] <- list(logl=out[[i]]$logl,coef=out[[i]]$coef,tau=out[[i]]$tau,pi=out[[i]]$pi,covar=out[[i]]$covar)
    }
  return(list(aic=aic,bic=bic,fm=fm))
 
}


"create.starting.values" <-
function (S,G,n,form,datsp) 
{
  ##apply.glm <- function(i,form,datsp,tau,n,dat.tau){
  ##  dat.tau <- rep(tau[,i],each=n)
  ##  dat.tau <- get("dat.tau")
  ##  datsp <- get("datsp")
    ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  ##  f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
    ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
   ## list(coef=f.mix$coef)
  ##}
  
  environment(form) <- environment()
  tau <- matrix(runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
  ##dat.tau <- rep(0,dim(datsp)[1])
 
  {fmM <- lapply(1:G,apply.glm,form,datsp,tau,n)}
  first.fit <- list(x=model.matrix(as.formula(form),data=datsp),y=datsp$obs)
 
  ##return(list(pi=pi,fmM=fmM,tau=tau,dat.tau=dat.tau,first.fit=first.fit))
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}


"distr.binom" <-
function( p){
 nobs <- length( p)
  new.dist <- old.dist <- rep( 0, nobs+1)
  old.dist[1] <- 1-p[1]
  old.dist[2] <- p[1]
  for( ii in 2:nobs){
    new.dist[1] <- old.dist[1]*(1-p[ii])
    for( jj in 2:ii)
      new.dist[jj] <- old.dist[jj-1]*p[ii] + old.dist[jj]*(1-p[ii])
    new.dist[ii+1] <- old.dist[ii]*p[ii]
    old.dist <- new.dist
  }
  return( new.dist)
}


"estimate.pi" <-
function (j,sp,spname,datsp,fmM,pi,G,first.fit) 
{
    link.fun <- make.link("logit")
    tmp.like <- rep(0,G)
  tau <- rep(0,G)
  sel.sp <- which(sp==spname[j])
  for(i in 1:G) {
    if(length(fmM[[i]]$coef)==1){lpre <- link.fun$linkinv(first.fit$x[sel.sp,]*fmM[[i]]$coef)
    }else{    lpre <- link.fun$linkinv(first.fit$x[sel.sp,]%*%fmM[[i]]$coef)}
   ##lpre <- link.fun$linkinv(first.fit$x[sel.sp,]%*%fmM[[i]]$coef)
    ## tmp.like[i] <- sum(dbinom(datsp$obs[sel.sp],1,fmM[[i]]$fitted[sel.sp],log=T))
    obs <- datsp$obs[sel.sp]
    lpre[obs==0]<- 1- lpre[obs==0]
    tmp.like[i] <- sum(log(lpre))
  }
  eps <- max(tmp.like)
  sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  for(i in 1:G) {
    tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
  }
  return(list(tau=tau,sum.like=sum.like))
}


"fitMix" <-
function (form,datsp,sp,G=2,ite.max=500,trace=TRUE,full.model=FALSE) 
{
## dat2 has colums obs,sp
  ##
  temp.warn <- getOption( "warn")
  options( warn=-1)

  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  cat("Fitting Group",G,"\n")
  if(trace) cat("Iteration | LogL \n")
  
  ##dat.tau <- data.frame(matrix(0,dim(datsp)[1],G))
  dat.tau <- 0
  ##dat <- data.frame(datsp,dat.tau)
  pi <- rep(0,G)
  ite <- 1
  logL <- -99999999
  old.logL <- -88888888

  ## set up initial GLM
  t1 <- create.starting.values(S,G,n,form,datsp)
  pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
 
 
  while(abs(logL-old.logL) > 0.0001 & ite<=ite.max){
    old.logL <- logL
    
    for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
    
    if(any(pi==0)) { ## occasionally with complicated models the random starting values result in a pi[i]==0; so restart with new random starts
      cat("pi has gone to zero - restarting fitting \n")
      t1 <- create.starting.values(S,G,n,form,datsp)
      pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
      ite <- 1
    }

 
    { fmM <- lapply(1:G,weighted.glm,first.fit,tau,n,fmM) }


    logL <- 0 
    tmp.like <- matrix(0,S,G)

    {est.tau <- lapply(1:S,estimate.pi,sp,sp.name,datsp,fmM,pi,G,first.fit)}
           
    for(j in 1:S){
      if(is.atomic(est.tau[[j]])){ print (est.tau[[j]])} else
      {
        tau[j,] <- est.tau[[j]]$tau
        logL <- logL+est.tau[[j]]$sum.like
      }
    }

    if(trace) cat(ite,"     | ",logL,"\n")
     ite <- ite+1
  }
  fm.out <- data.frame(matrix(0,G,length(fmM[[1]]$coef)))
  names(fm.out) <- names(fmM[[1]]$coef)
  tau <- data.frame(tau)
  names(tau) <- paste("grp.",1:G,sep="")
  EN <- -sum(unlist(tau)*log(unlist(tau)))
  d <- length(unlist(fm.out)) + length(tau)-1
  for(i in 1:G) {
    fm.out[i,] <- fmM[[i]]$coef
    ##dat.tau[,i] <- rep(tau[,i],each=n)
  }

  names(pi) <- paste("G",1:G,sep=".")
  t.pi <- additive.logistic(pi,TRUE)
  parms <- c(t.pi[1:(G-1)],unlist(fm.out))
  logL.full <- logL
  logL <- logLmix(parms,first.fit,G,S,sp,sp.name)

  options(warn=temp.warn)
  if(full.model)  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,fmM=fmM,model.tau=dat.tau,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
  
  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
  
}


"fitmix.cpp" <-
function (form,datsp,sp,G=2,pars=NA,trace=TRUE,calc.hes=FALSE) 
{
  if(!is.numeric(sp)){
    sp <- as.integer(factor(sp))
  }
  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  X <- model.matrix(form, data = datsp[sp==sp.name[1],])
  ##X <- model.frame(form, data = datsp)
 # X <- X[sp==sp.name[1],]
  
  #y <- model.response(form, data = datsp)
  y <- datsp$obs
  if(is.na(pars[1])) {
    ##pars <- rep(0.01,G-1+(ncol(X)*G))
    pars <- runif(G-1+(ncol(X)*G),-2,2)
    pars[1:(G-1)] <- runif(G-1)
  }
##  hes <- rep(0,length(pars)^2)
  gradient <- rep(0,length(pars))
  tau <- matrix(0,S,G) ##must leave this in as defines S & G
  
  loglike <- try(.Call("SpeciesMix",pars,y,X,sp,tau,gradient,PACKAGE="SpeciesMix"))
  
  calc.deriv <- function(p){
    gradient <- rep(0,length(pars))
    ll <- .Call("Calculate_Gradient",p,y,X,sp,tau,gradient,PACKAGE="SpeciesMix")
    return(gradient)
  }
  hes <- 0
  if(calc.hes){
    hes <- nd2(pars,calc.deriv)
    dim(hes) <- rep(length(pars),2)
    dim(hes) <- rep(length(pars),2)
    rownames(hes) <- colnames(hes) <- c(paste("G.",1:(G-1),sep=""),paste("G",1:G,rep(colnames(X),each=G),sep="."))
  }
  if(!is.numeric(loglike)) loglike <- 0
  pi <- pars[1:(G-1)]
  coef <- pars[ (G):length(pars)]

  r.logl <- logLmix(pars,list(y=y,x=model.matrix(form, data = datsp)),G,S,sp,sp.name,out.tau=TRUE)
  pi <- additive.logistic(pi)
  names(pi) <- paste("G.",1:G,sep="")
  coef <- matrix(coef,G,ncol(X))
  rownames(coef) <- paste("G.",1:G,sep="")
  colnames(coef) <- colnames(X)

  AIC <- 2*loglike + 2*length(pars)
  BIC <- 2*loglike + log(S)*length(pars)
  ##list(logl=loglike,r.logl=r.logl$logl,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient)
  list(logl=loglike,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient)
  
}


"Lmix" <-
function (pars,y,x,G) 
{
   fm <- pars[-1*(1:(G-1))]
    pi <- pars[(1:(G-1))]
    dim(fm) <- c(G,(length(pars)-(G-1))/G)
    pi <- additive.logistic(pi)
    S <- dim(y)[2]
   
   link.fun <- make.link("logit")

    like <- 1
    lf <- matrix(0,G,S)
    lg <- rep(0,S)
    dBi <- array(0,dim=c(S,G,dim(fm)[2]))
    for(s in 1:S){
      tmp.like <- rep(0,G)
      for(g in 1:G){
        lpre <- x%*%fm[g,]
        for(j in 1:dim(fm)[2]){
          dBi[s,g,j] <- sum((y[,s]-link.fun$linkinv(lpre))*x[,j])
        }
        tmp.like[g] <- pi[g]*prod(dbinom(y[,s],1,link.fun$linkinv(lpre),log=F))
        lf[g,s] <- prod(dbinom(y[,s],1,link.fun$linkinv(lpre),log=F))
      }
      lg[s] <- sum(tmp.like)
      like <- like*sum(tmp.like)
    }
    #print(dBi)
  #  print(log(lg))
  #  print(log(lf))
dl.dpi <- rep(0,G)
   for(g in 1:G) dl.dpi[g] <- ( sum( exp(-log(lg) + log(lf[g,]))))
   
    der <- matrix(0,dim(fm)[1],dim(fm)[2])
    for(j in 1:dim(fm)[2]){
      for(g in 1:G){
        for(s in 1:S){
          der[g,j] <- der[g,j]+ exp( -log(lg[s]) + log(pi[g]) + log(lf[g,s])) * dBi[s,g,j]
          #if(g==1 & j == 1) print(c(der[g,j],-log(lg[s]),pi[g],log(lf[g,s]),dBi[s,g,j]))
        }
      }
    }
dpi.deta <- matrix(0,G-1,G)
   ad.trans <- 1+sum(exp(pars[1:(G-1)]))
   for(i in 1:(G-1))
     for(g in 1:G-1){
       if(i==g) {dpi.deta[i,g] <- exp(pars[i])*exp(pars[g])/ad.trans^2}
       else{ dpi.deta[i,g] <- exp(pars[i])/ad.trans - exp(2*pars[g])/ad.trans^2}
     }
   dpi.deta[,G] <- -rowSums(dpi.deta)
   print(dpi.deta)
    print(dl.dpi)
   d1 <- dpi.deta%*%dl.dpi
   print(d1)
    
    list(like,0-der)
}


"logLike.pars" <-
function (pi,coef,sp.form,sp.data,covar.data) 
{

  G <- length(pi)
  pars <- c(additive.logistic(pi,T)[1:(G-1)],unlist(coef))
  
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data))
  logl <- -logLmix(pars,first.fit,G,S,sp,sp.name)
  logl
}


"logLmix" <-
function (pars,first.fit,G,S,sp,spname,out.tau=FALSE) 
{
   tau <- matrix(0,S,G)
##tau,out.tau=FALSE
  if(G>1){
    fm <- pars[-1*(1:(G-1))]
    pi <- pars[(1:(G-1))]
##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
    #pi <- tau[((length(tau)-(G-2)):length(tau))]
    dim(fm) <- c(G,(length(pars)-(G-1))/G)
    ##pi[G] <- 1-sum(pi)
    pi <- additive.logistic(pi)
  } else{
    fm <- tau[1:(length(pars)-1)]
    dim(fm) <- c(1,length(fm))
    pi <- 1
  }
  
  link.fun <- make.link("logit")

  
 log.like <- 0
  for(j in 1:S){
    sel.sp <- which(sp==spname[j])
    tmp.like <- rep(0,G)
    for(i in 1:G){
      if(length(fm[i,])==1){lpre <- first.fit$x[sel.sp,]*fm[i,]
      }else{      lpre <- first.fit$x[sel.sp,]%*%fm[i,]}
      tmp.like[i] <- sum(dbinom(first.fit$y[sel.sp],1,link.fun$linkinv(lpre),log=T))
    }
    eps <- max(tmp.like)
    log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
    tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  }
  
  if(out.tau)return(list(logl=log.like,tau=tau))
  log.like
}


"mix.residuals" <-
function (fmM,form,datsp,sp) 
{
   cat("calculating residuals \n")
  link.fun <- make.link("logit")
  x <- model.matrix(form,data=datsp)
  spname <- unique(sp)
  S <- length(spname)
  G <- ncol(fmM$tau)

  PIT <- matrix(NA,S,G)
  for(g in 1:G){
    for(s in 1:length(spname)){
      sel.sp <- which(sp==spname[s])
      t.obs <- sum(datsp$obs[sel.sp])
      pre <- link.fun$linkinv(x[sel.sp,]%*%fmM$coef[g,])
      obs <- datsp$obs[sel.sp]
      dis <- distr.binom(pre)
#      PIT[s,g] <- qnorm(sum(dis[1:(length(dis)-1)],dis[length(dis)]/2))
      nSucc <- sum( obs)
      transfo <- sum( dis[1:nSucc],dis[nSucc+1]/2)
      transfo <- min( transfo, 1)
      ##transfor <- max( transfo, 0)
      PIT[s,g] <- qnorm( transfo)
    }
  }
  PIT 
}


"nd2" <-
function( x0, f, m=NULL, D.accur=4, ...) {
# A function to compute highly accurate first-order derivatives 
# From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)  
# Adapted by Scott Foster from code nicked off the net 2007

  D.n <- length( x0)
  if ( is.null( m)) {
    D.f0 <- f(x0, ...)
    m <- length( D.f0) 
  }
  if ( D.accur == 2) {
    D.w <- tcrossprod( rep( 1, m),c( -1/2, 1/2))
    D.co <- c( -1, 1) 
  }
  else {
    D.w <- tcrossprod( rep( 1, m),c( 1/12, -2/3, 2/3, -1/12))
    D.co <- c( -2, -1, 1, 2) 
  }
  D.n.c <- length( D.co)
  macheps <- .Machine$double.eps
  D.h <- macheps^( 1/3)*abs( x0)
  D.deriv <- matrix( NA, nrow=m, ncol=D.n)
  for ( ii in 1:D.n) {
    D.temp.f <- matrix( 0, m, D.n.c)
    for ( jj in 1:D.n.c) {
      D.xd <- x0+D.h[ii]*D.co[jj]*( 1:D.n == ii)
      D.temp.f[,jj] <- f( D.xd, ...) 
    }
    D.deriv[,ii] <- rowSums( D.w*D.temp.f)/D.h[ii] 
  }
  return( as.double( D.deriv))
}


"nH2" <-
function( pt, fun, accur=c(4,4), type="H.Diag", ...) {
# A function to compute highly accurate second order Hessian diags and other off-diags
# partially from Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
# Adapted by Scott Foster from code nicked off the net 2007

  H.n <- length( pt)
  derivs <- function( d.x0, ...) { nd2( x0=d.x0, f=fun, m=1, D.accur=accur[2], ...) }
  Hes <- nd2( x0=pt, f=derivs, D.accur=accur[2], ...)
  Hes <- matrix( Hes, nrow=length( pt))
  Hes <- ( Hes+t( Hes))/2

  if ( type == "H.Diag") {
    macheps <- .Machine$double.eps
    H.h <- macheps^(1/4)*abs(pt)
    H.f0 <- fun( pt, ...)
    H.m <- length( H.f0)
    if ( accur[1] == 2) {
      H.w <- tcrossprod( rep( 1, H.m), c( 1, -2, 1))
      H.co <- c( -1, 0, 1) 
    }
    else {
      H.w <- tcrossprod( rep( 1, H.m), c( -1/12, 4/3, -5/2, 4/3, -1/12))
      H.co <- c( -2, -1, 0, 1, 2) 
    }
    H.n.c <- length( H.co)
    Hes.diag <- double( length=H.n)
    for ( ii in 1:H.n) {
      H.temp.f <- matrix( 0, H.m, H.n.c)
      for ( jj in 1:H.n.c) {
        if ( H.co[jj] != 0) {
          H.xd <- pt+H.h[ii]*H.co[jj]*( 1:H.n == ii)
          H.temp.f[,jj] <- fun( H.xd, ...) 
        }
        else
          H.temp.f[,jj] <- H.f0 
      }
      Hes.diag[ii] <- rowSums( H.w*H.temp.f)/( H.h[ii]^2) 
    } 
    diag( Hes) <- Hes.diag 
  }
  return( Hes)
}


".onLoad" <-
function (libname, pkgname) 
{
  require(MASS)
     # Generic DLL loader
    dll.path <- file.path( libname, pkgname, 'libs')
    if( nzchar( subarch <- .Platform$r_arch))
      dll.path <- file.path( dll.path, subarch)
    this.ext <- paste( sub( '.', '[.]', .Platform$dynlib.ext, 
fixed=TRUE), '$', sep='')

    dlls <- dir( dll.path, pattern=this.ext, full=FALSE)
    names( dlls) <- dlls
    if( length( dlls))
      lapply( dlls, function( x) library.dynam( sub( this.ext, '', x), 
package=pkgname, lib.loc=libname))

}


"plot.archetype" <-
function (mixture.model,pch=20,...) 
{
  matplot(mixture.model$residuals*mixture.model$tau,pch=pch,col=1:ncol(mixture.model$tau),...,xlab="Species",ylab="Residuals * tau")
}


"predict.archetype" <-
function (object,new.obs,...)
{
  mixture.model <- object
  G <- length(mixture.model$pi)
  covar <- mixture.model$covar[-(1:(G-1)),-(1:(G-1))]
  coef <- mixture.model$coef
  model.fm <- as.formula(mixture.model$formula)
  model.fm[[2]] <- NULL
  X <- model.matrix(model.fm,new.obs)
  link.fun <- make.link("logit")
  outvar <- matrix(NA,dim(X)[1],G)
  outpred <- matrix(NA,dim(X)[1],G)
  colnames(outvar) <- colnames(outpred) <- paste("G",1:G,sep=".")
  for(g in 1:G){
    lp <- as.numeric(X%*%coef[g,])
    outpred[,g] <- link.fun$linkinv(lp)
    dhdB <- (exp(lp)/(1+exp(lp)))*X - exp(lp)^2/((1+exp(lp))^2)*X
    c2 <- covar[seq(g,dim(covar)[1],G),seq(g,dim(covar)[1],G)]
    ##pos <- seq(1,dim(X)[1],1000)
    ##for(k in pos[-length(pos)]){
    ##  outvar[k:(k+1000),g] <- diag(tcrossprod(dhdB[k:(k+1000),]%*%c2,dhdB[k:(k+1000),]))
    ##}
    ##outvar[pos[length(pos)]:dim(X)[1],] <- diag(tcrossprod(dhdB[pos[length(pos)]:dim(X)[1],]%*%c2,dhdB[pos[length(pos)]:dim(X)[1],]))
    ##outvar[,g] <- diag(tcrossprod(dhdB%*%c2,dhdB))
    for(k in 1:dim(X)[1]){
      outvar[k,g] <- (dhdB[k,]%*%c2)%*%(dhdB[k,])
    }
  }
  list(fit=outpred,se.fit=sqrt(outvar))
}




"SpeciesMix" <-
function (sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=4, em.refit = 1 , est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  if(em.prefit | G==1){
    prefit <- SpeciesMix.em(sp.form,sp.data,covar.data,G,ite.max=em.steps,em.refit=em.refit)
    pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
  }
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  if(G==1) return(prefit)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  while(fmM.out$logl==0) {
    prefit <- SpeciesMix.em(sp.form,t.sp.data,t.covar.data,G,ite.max=em.steps,em.refit=1)
    pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
    fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  }

  rownames(fmM.out$tau) <- sp.name ## add names to taus

  if(est.var){
    fmM.out$covar <- try(solve(fmM.out$hessian))
  }
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form
  class(fmM.out) <- "archetype"
  fmM.out
}



"SpeciesMix.em" <-
function (sp.form,sp.data,covar.data,G=2,em.refit=1,ite.max=500, est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  
  fmM.out <- fitMix(sp.form,data,sp,G,ite.max,trace)
  if(em.refit>1)
    for(i in 2:em.refit){
      fmM <- fitMix(sp.form,data,sp,G,ite.max,trace)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  
  if(est.var){
    var <- 0
    t.pi <- additive.logistic(fmM.out$pi,TRUE)
    parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef))
    first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data))
    fun.est.var <- function(x){-logLmix(x,first.fit,G,S,sp,sp.name)}
    var <- solve( nH2( pt=parms, fun=fun.est.var))
    colnames( var) <- rownames( var) <- names( parms)
    fmM.out$covar <- var
  }
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}




"weighted.glm" <-
function (g,first.fit,tau,n,fmM) 
{
  dat.tau <- rep(tau[,g],each=n)
  ##lpre <- first.fit$x%*%fmM[[g]]$coef
  ##f.mix <- glm.fit(x=first.fit$x,y=first.fit$y,weights=dat.tau[,g],family=binomial(),etastart=fmM[[g]]$linear.predictors)
  f.mix <- glm.fit(x=first.fit$x,y=first.fit$y,weights=dat.tau,family=binomial(),start=fmM[[g]]$coef)
  ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
  ##list(linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
  list(coef=f.mix$coef)
  
}

