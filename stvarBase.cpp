#include"stvarBase.hpp"
#include<cmath>
#include<R.h>
#include<Rmath.h>
stvarBase::stvarBase(Mmat xmat_, Mmat resp_, int p_, const VectorXd& delta,
		     double s0_, double n0_, int nmcmc_, int nburn_,
		     int nthin_, int nrec_, double powrho, double powres,
		     int priRho, const VectorXd& priParRho,
		     int priRes, const VectorXd& priParRes,
		     int kerRho, const VectorXd& parKerRho,
		     int kerRes, const VectorXd& parKerRes,
		     const VectorXd& thetaRes0):
  xmat(xmat_.data(),xmat_.rows(),xmat_.cols()), resp(resp_.data(),resp_.rows(),resp_.cols()),
  p(p_), s0(s0_), n0(n0_), nmcmc(nmcmc_), nburn(nburn_), nthin(nthin_), nrec(nrec_)
{
  n = xmat_.rows();
  d = xmat_.cols();
  T = resp_.cols();
  Tinf = T - p;
  dobs = delta(0);
  dsys = delta(1);
  vVar.resize(Tinf), st.resize(Tinf);
  nt.resize(Tinf), Err.resize(Tinf);
  spresp = VectorXd::Zero(n);
  for(int i = p; i != T; ++i)
    spresp += resp.col(i);
  spresp /= Tinf;
  
  pcorrRho = new powerexpCorr(d,VectorXd::Constant(d,powrho));
  pcorrRes = new powerexpCorr(d,VectorXd::Constant(d,powres));
  setPriorRho(priRho,priParRho);
  setPriorRes(priRes,priParRes);
  setKernelRho(kerRho,parKerRho);
  setKernelRes(kerRes,parKerRes);
  initThetaRes(thetaRes0);
}

stvarBase::stvarBase(Mmat xmat_, Mmat resp_, int p_, const VectorXd& delta,
		     double s0_, double n0_, int nmcmc_, int nburn_,
		     int nthin_, int nrec_, double powres,
		     int priRes, const VectorXd& priParRes,
		     int kerRes, const VectorXd& parKerRes,
		     const VectorXd& thetaRes0):
  xmat(xmat_.data(),xmat_.rows(),xmat_.cols()), resp(resp_.data(),resp_.rows(),resp_.cols()),
  p(p_), s0(s0_), n0(n0_), nmcmc(nmcmc_), nburn(nburn_), nthin(nthin_), nrec(nrec_)
{
  n = xmat_.rows();
  d = xmat_.cols();
  T = resp_.cols();
  Tinf = T - p;
  dobs = delta(0);
  dsys = delta(1);
  vVar.resize(Tinf), st.resize(Tinf);
  nt.resize(Tinf), Err.resize(Tinf);
  pcorrRes = new powerexpCorr(d,VectorXd::Constant(d,powres));
  setPriorRes(priRes,priParRes);
  setKernelRes(kerRes,parKerRes);
  initThetaRes(thetaRes0);
}

void stvarBase::setPriorRho(int priRho, const VectorXd& priParRho)

{
  switch(priRho)
  {
  case GAMMA:
    priThetaRho = new gammaDistb(priParRho(0), priParRho(1));
    break;
  case INVGA:
    priThetaRho = new invgammaDistb(priParRho(0), priParRho(1));
    break;
  default:
    error("unrecognized prior type for Rho correlation!");
  }
}
void stvarBase::setPriorRes(int priRes, const VectorXd& priParRes)

{
  switch(priRes)
  {
  case GAMMA:
    priThetaRes = new gammaDistb(priParRes(0), priParRes(1));
    break;
  case INVGA:
    priThetaRes = new invgammaDistb(priParRes(0), priParRes(1));
    break;
  default:
    error("unrecognized prior type for Rho correlation!");
  }
}

void stvarBase::setKernelRho(int kerRho, const VectorXd& parKerRho)
			     
{
   switch(kerRho)
   {
   case UNIF:
     pKerThetaRho = new uniformKernel(parKerRho(0),parKerRho(1),
				      parKerRho(2),parKerRho(3)>0.0);
     break;
   case LNORMAL:
     pKerThetaRho = new lognormalKernel(parKerRho(0));
    break;
   case NORMAL:
     pKerThetaRho = new normalKernel(parKerRho(0));
   default:
     error("unrecognized kernel type for correlation of rho");
   }
}

void stvarBase::setKernelRes(int kerRes, const VectorXd& parKerRes)
			     
{
   switch(kerRes)
   {
   case UNIF:
     pKerThetaRes = new uniformKernel(parKerRes(0),parKerRes(1),
				      parKerRes(2),parKerRes(3)>0.0);
     break;
   case LNORMAL:
     pKerThetaRes = new lognormalKernel(parKerRes(0));
    break;
   case NORMAL:
     pKerThetaRes = new normalKernel(parKerRes(0));
   default:
     error("unrecognized kernel type for correlation of residual");
   }
}

void stvarBase::initThetaRes(const VectorXd& theta0)
{
  thetaRes = theta0;
  MatrixXd corr0(n,n);
  pcorrRes -> evalCorrMat(xmat,theta0,corr0);
  lltcorrRes.compute(corr0);
  if(lltcorrRes.info() != Eigen::Success)
    error("Cholesky failure while initializing correlation of residuals!");
  corrRes = corr0;
}
void stvarBase::evalLoglikRes()
{
  MatrixXd lmat = lltcorrRes.matrixL();
  double logdet = 2.0*lmat.diagonal().array().log().sum();
  loglikRes = -1.0*Tinf*logdet;
  VectorXd err;
  double vt;
  for(int i = 0; i != Tinf; ++i)
  {
    vt = vVar[i];
    err = Err[i];
    loglikRes -= err.dot(lltcorrRes.solve(err))/vt;
  }
}

double stvarBase::evalNewLoglikRes(const VectorXd& newtheta, MatrixXd& newcorr,
				   LLTMd& newllt)
{
  newcorr = MatrixXd::Zero(n,n);
  pcorrRes -> evalCorrMat(xmat, newtheta, newcorr);
  newllt.compute(newcorr);
  if(newllt.info() != Eigen::Success)
    error("Cholesky failure while evaluating new log likelihood for residuals!");
  MatrixXd lmat = newllt.matrixL();
  double logdet = 2.0*lmat.diagonal().array().log().sum();
  double loglik = -1.0*Tinf*logdet;
  VectorXd err;
  double vt;
  for(int i = 0; i != Tinf; ++i)
  {
    vt = vVar[i];
    err = Err[i];
    loglik -= err.dot(newllt.solve(err))/vt;
  }
  return loglik;
}

void stvarBase::drawThetaRes()
{
  VectorXd newtheta(d);
  double newthetai, newlogpri = 0.0;
  double logforden = 0.0, logbackden = 0.0;
  for(int i = 0; i != d; ++i)
  {
    newthetai = pKerThetaRes -> draw(thetaRes(i));
    newlogpri += priThetaRes -> logpdf(newthetai);
    logforden += pKerThetaRes -> logdensity(thetaRes(i),newthetai);
    logbackden += pKerThetaRes -> logdensity(newthetai,thetaRes(i));
    newtheta(i) = newthetai;
  }
  MatrixXd newcorr;
  LLTMd newllt;
  double newloglik = evalNewLoglikRes(newtheta, newcorr, newllt);
  double logpost = loglikRes + logPriRes;
  double newlogpost = newloglik + newlogpri;
  double logratio = newlogpost - logpost + logbackden - logforden;
  double u = runif(0.0,1.0);
  if(log(u)<logratio)
  {
    thetaRes = newtheta;
    lltcorrRes = newllt;
    corrRes = newcorr;
    loglikRes = newloglik;
    logPriRes = newlogpri;
  }
}
double stvarBase::evalLogLikRho(const VectorXd& theta, const VectorXd& rho,
				MatrixXd& newcorr, LLTMd& lltofCorr)
{
  MatrixXd lmat;
  VectorXd solrho;
  double loglik, quadrho;
  pcorrRho -> evalCorrMat(xmat,theta,newcorr);
  lltofCorr.compute(newcorr);
  if(lltofCorr.info() != Eigen::Success)
    error("Cholesky failure while evalulating likelihood function for rho!");
  solrho = lltofCorr.solve(rho);
  quadrho = rho.dot(solrho);
  lmat = lltofCorr.matrixL();
  loglik = -lmat.diagonal().array().log().sum() - 0.5*quadrho;
  return loglik;
}

void stvarBase::drawVar()
{
  int idx = Tinf-1;
  double frac = 1.0 - dobs;
  double niter = nt[idx];
  double diter = niter*st[idx];
  double viter = rgamma(0.5*niter,2.0/diter);
  vVar[idx] = 1.0/viter;
  for(--idx ;idx >= 0; --idx)
  {
    niter = frac * nt[idx];
    diter = niter * st[idx];
    viter = dobs * viter + rgamma(0.5*niter,2.0/diter);
    vVar[idx] = 1.0/viter;
  }
}

void stvarBase::mvNormalDrawPrec(const VectorXd& likmean, const MatrixXd& prec,
				 double sigma, int dim, VectorXd& out)
{
  LLTMd lltofprec(prec);
  if(lltofprec.info() != Eigen::Success)
    error("Cholesky failure for precision matrix in multivariate normal sampling!");

  for(int i = 0; i != dim; ++i)
    out(i) = norm_rand();
  out.noalias() = lltofprec.matrixU().solve(out)*sigma;
  out += lltofprec.solve(likmean);
}
