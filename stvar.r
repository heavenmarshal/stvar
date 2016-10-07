#dyn.load("stvar.so")
stvar <- function(xx,resp,p,newxx,newystart,
                  delta,s0=10,n0=1,nmcmc=5500,nburn=500,nthin=10,
                  powrho=2,powres=1.999,priRho=c("gamma","invgamma"),
                  priRes=c("gamma","invgamma"), kerRho=c("lnormal","normal","unif"),
                  parKerRho, kerRes=c("lnormal","normal","unif"), parKerRes,
                  C0=10*diag(p), m0 = rep(1,p))
{
    nn <- nrow(xx)
    if(nrow(resp) != nn) stop("the dimension of input and response matrix do not meet!")
    dd <- ncol(xx)
    if(ncol(newxx) != dd) stop("the dimension of training and prediction input matrices do not meet!")
    if(ncol(newystart) != p) stop("the starting value of prediction responses must have p time stages!")
    TT <- ncol(resp)
    if(length(delta) < 2 ) stop("the length of distcounting parameter vector can not be smaller than 2!")
    if(nmcmc<nburn) stop("nmcmc must be no less than nburn!")
    if((nmcmc-nburn)%%nthin != 0)
    {
        nmcmc <- nmcmc - (nmcmc-nburn)%%nthin
        warning("the number of mcmc iterations has been reduced to ",nmcmc)
    }
    nrec <- (nmcmc-nburn)/nthin+1
    newn <- nrow(newxx)
    outlen <- nrec*newn*TT
    priRho <- match.arg(priRho)
    priRho <- switch(priRho,
                     gamma = 101,
                     invgamma = 102)
    priRes <- match.arg(priRes)
    priRes <- switch(priRes,
                     gamma = 101,
                     invgamma = 102)
    priorinfo <- laGP::darg(NULL,xx)
    priParRes <- priParRho <- priorinfo$ab
    if(priRho == 102) priParRho[2] <- 1/priParRho[2]
    if(priRes == 102) priParRes[2] <- 1/priParRes[2]
    thetaRho0 <- thetaRes0 <- rep(priorinfo$start,dd)
    kerRho <- match.arg(kerRho)
    lenkerrho <- switch(kerRho,
                        lnormal = 1,
                        normal = 1,
                        unif = 4)
    if(kerRho=="unif")
        parKerRho <- c(parKerRho,min(.5,priorinfo$max),priorinfo$min,FALSE)
    if(length(parKerRho) != lenkerrho)
        stop("incorrect length of kernel parameters for Rho!")
    kerRho <- switch(kerRho,
                     lnormal = 202,
                     normal = 203,
                     unif = 201)
    kerRes <- match.arg(kerRes)
    lenkerres <- switch(kerRes,
                        lnormal = 1,
                        normal = 1,
                        unif = 4)
    if(kerRes=="unif")
        parKerRes <- c(parKerRes,priorinfo$max,priorinfo$min,FALSE)
    if(length(parKerRes) != lenkerres)
        stop("incorrect length of kernel parameters for Residuals!")
    kerRes <- switch(kerRes,
                     lnormal = 202,
                     normal = 203,
                     unif = 201)
    set.seed(as.integer(Sys.time()))
    out <- .C("stvar_R",as.double(xx),as.double(resp),as.integer(nn),as.integer(TT),
              as.integer(dd),as.integer(p),as.double(delta),as.double(s0),as.double(n0),
              as.integer(nmcmc),as.integer(nburn),as.integer(nthin),as.integer(nrec),
              as.double(powrho),as.double(powres),as.integer(priRho),as.double(priParRho),
              as.integer(2),as.integer(priRes),as.double(priParRes),as.integer(2),
              as.integer(kerRho), as.double(parKerRho),as.integer(lenkerrho),
              as.integer(kerRes), as.double(parKerRes),as.integer(lenkerres),
              as.double(thetaRho0), as.double(thetaRes0),
              as.double(C0),as.double(m0),as.double(newxx),as.double(newystart),
              as.integer(newn), pred=double(outlen),rho=double(nrec*p*newn),
              thetaRho=double(nrec*p*dd), phi = double((TT-p)*p*nrec))
    pred <- array(out$pred,c(newn,TT,nrec))
    rho <- array(out$rho,c(newn,p,nrec))
    thetaRho <- array(out$thetaRho,c(dd,p,nrec))
    phi <- array(out$phi,c(p,TT-p,nrec))
    return(list(pred=pred,rho=rho,thetaRho=thetaRho, phi=phi))
}

tvar <- function(xx,resp,p,newxx,newystart,
                 delta,s0=10,n0=1,nmcmc=5500,nburn=500,nthin=10,
                 powres=1.999, priRes=c("gamma","invgamma"), kerRes=c("lnormal","normal","unif"),
                 parKerRes, C0=10*diag(p), m0 = rep(1,p))
{
    nn <- nrow(xx)
    if(nrow(resp) != nn) stop("the dimension of input and response matrix do not meet!")
    dd <- ncol(xx)
    if(ncol(newxx) != dd) stop("the dimension of training and prediction input matrices do not meet!")
    if(ncol(newystart) != p) stop("the starting value of prediction responses must have p time stages!")
    TT <- ncol(resp)
    if(length(delta) < 2 ) stop("the length of distcounting parameter vector can not be smaller than 2!")
    if(nmcmc<nburn) stop("nmcmc must be no less than nburn!")
    if((nmcmc-nburn)%%nthin != 0)
    {
        nmcmc <- nmcmc - (nmcmc-nburn)%%nthin
        warning("the number of mcmc iterations has been reduced to ",nmcmc)
    }
    nrec <- (nmcmc-nburn)/nthin+1
    newn <- nrow(newxx)
    outlen <- nrec*newn*TT
    priRes <- match.arg(priRes)
    priRes <- switch(priRes,
                     gamma = 101,
                     invgamma = 102)
    priorinfo <- laGP::darg(NULL,xx)
    priParRes <- priorinfo$ab
    if(priRes == 102) priParRes[2] <- 1/priParRes[2]
    kerRes <- match.arg(kerRes)
    lenkerres <- switch(kerRes,
                        lnormal = 1,
                        normal = 1,
                        unif = 4)
    thetaRes0 <- rep(priorinfo$start,dd)
    if(length(parKerRes) != lenkerres)
        stop("incorrect length of kernel parameters for Residuals!")
    kerRes <- switch(kerRes,
                     lnormal = 202,
                     normal = 203,
                     unif = 201)
    set.seed(as.integer(Sys.time()))
    out <- .C("tvar_R",as.double(xx),as.double(resp),as.integer(nn),as.integer(TT),
              as.integer(dd),as.integer(p),as.double(delta),as.double(s0),as.double(n0),
              as.integer(nmcmc),as.integer(nburn),as.integer(nthin),as.integer(nrec),
              as.double(powres), as.integer(priRes),as.double(priParRes),as.integer(2),
              as.integer(kerRes), as.double(parKerRes),as.integer(lenkerres),
              as.double(thetaRes0), as.double(C0),as.double(m0),as.double(newxx),
              as.double(newystart), as.integer(newn), pred=double(outlen))
    pred <- array(out$pred,c(newn,TT,nrec))
    return(pred)
}
