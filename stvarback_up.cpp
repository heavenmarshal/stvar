#include"stvar.hpp"
#include<cmath>
#include<R.h>
#include<Rmath.h>
#include<algorithm>

using Eigen::TriangularView;
using Eigen::Upper;
using Eigen::Lower;
typedef TriangularView<MatrixXd, Upper> ViewUmd;
typedef TriangularView<const MatrixXd, Lower> ViewLcmd;

void stvar::initThetaRho(const VectorXd& theta0)
{
  MatrixXd corrrho(n,n);
  pcorrRho -> evalCorrMat(xmat,theta0,corrrho);
  LLTMd lltcorrrho(corrrho);
  if(lltcorrrho.info() != Eigen::Success)
    error("Cholesky failure while initializing correlation for rho!");
  for(int i = 0; i != p; ++i)
  {
    thetaRho.push_back(theta0);
    lltRho.push_back(lltcorrrho);
  }
}
void stvar::forwardFilter()
{
  MatrixXd F, R, A, Ad, Q, C=C0;
  VectorXd yvec, err, mvec = m0;
  double siter = s0, niter = n0, diter;
  LLTMd lltofq;
  for(int i = p, j=0; i != T; ++i,++j)
  {
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    Ad = R*F.transpose();
    Q = F*Ad + siter*corrRes;
    lltofq.compute(Q);
    if(lltofq.info() != Eigen::Success)
      error("Cholesky failure while filtering forward!");
    A = lltofq.solve(Ad.transpose()).transpose();
    mvec += A*err;
    C = R-A*Q*A.transpose();
    C /= siter;
    diter = dobs * niter * siter;
    diter += siter*err.dot(lltofq.solve(err));
    niter = dobs * niter + n;
    siter = diter/niter;
    C *= siter;
    C.noalias() = 0.5*(C+C.transpose());
    Ct[j] = C;
    mt[j] = mvec;
    nt[j] = niter;
    st[j] = siter;
  }
}
void stvar::tvarForwardFilter()
{
  MatrixXd F, R, A, Ad, Q, C=C0;
  VectorXd yvec, err, mvec = m0;
  double siter = s0, niter = n0, diter;
  LLTMd lltofq;
  for(int i = p, j=0; i != T; ++i,++j)
  {
    F = resp.block(0,j,n,p);
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    Ad = R*F.transpose();
    Q = F*Ad + siter*corrRes;
    lltofq.compute(Q);
    if(lltofq.info() != Eigen::Success)
      error("Cholesky failure while filtering forward!");
    A = lltofq.solve(Ad.transpose()).transpose();
    mvec += A*err;
    C = R-A*Q*A.transpose();
    C /= siter;
    diter = dobs * niter * siter;
    diter += siter*err.dot(lltofq.solve(err));
    niter = dobs * niter + n;
    siter = diter/niter;
    C *= siter;
    C.noalias() = 0.5*(C+C.transpose());
    Ct[j] = C;
    mt[j] = mvec;
    nt[j] = niter;
    st[j] = siter;
  }
}
void stvar::forwardSmooth()
{
  MatrixXd F, R, A, Ad, Q, C = C0;
  VectorXd yvec, err, mvec = m0;
  LLTMd lltofq;
  for(int i = p, j = 0; i != T; ++i, ++j)
  {
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    Ad = R*F.transpose();
    Q = F*Ad + vVar[j]*corrRes;
    lltofq.compute(Q);
    if(lltofq.info() != Eigen::Success)
      error("Cholesky failure while smoothing forward!");
    A = lltofq.solve(Ad.transpose()).transpose();
    mvec += A*err;
    C = R - A*Q*A.transpose();
    C.noalias() = 0.5*(C+C.transpose());
    mt[j] = mvec;
    Ct[j] = C;
  }
}
void stvar::tvarForwardSmooth()
{
  MatrixXd F, R, A, Ad, Q, C = C0;
  VectorXd yvec, err, mvec = m0;
  LLTMd lltofq;
  for(int i = p, j = 0; i != T; ++i, ++j)
  {
    F = resp.block(0,j,n,p);
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    Ad = R*F.transpose();
    Q = F*Ad + vVar[j]*corrRes;
    lltofq.compute(Q);
    if(lltofq.info() != Eigen::Success)
      error("Cholesky failure while smoothing forward!");
    A = lltofq.solve(Ad.transpose()).transpose();
    mvec += A*err;
    C = R - A*Q*A.transpose();
    C.noalias() = 0.5*(C+C.transpose());
    mt[j] = mvec;
    Ct[j] = C;
  }
}

void stvar::backwardSample()
{
  int i = T-1, j = Tinf-1;
  VectorXd phi = VectorXd::Zero(p);
  VectorXd mvec = mt[j];
  MatrixXd C = Ct[j];
  VectorXd yvec = resp.col(i);
  MatrixXd F = Ft[j];
  VectorXd err;
  double frac = 1.0-dobs;
  mvNormalDrawCov(mvec, C, 1.0, p, phi);
  err = yvec - F*phi;
  Err[j] = err;
  vPhi[j] = phi;
  for(--i, --j; j >=0; --i, --j)
  {
    mvec = frac * mt[j] + dobs * phi;
    C = frac * Ct[j].array();
    mvNormalDrawCov(mvec,C,1.0,p,phi);
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*phi;
    Err[j] = err;
    vPhi[j] = phi;
  }
}

void stvar::tvarBackwardSample()
{
  int i = T-1, j = Tinf-1;
  VectorXd phi = VectorXd::Zero(p);
  VectorXd mvec = mt[j];
  MatrixXd C = Ct[j];
  VectorXd yvec = resp.col(i);
  MatrixXd F = resp.block(0,j,n,p);
  VectorXd err;
  double frac = 1.0-dobs;
  mvNormalDrawCov(mvec, C, 1.0, p, phi);
  err = yvec - F*phi;
  Err[j] = err;
  vPhi[j] = phi;
  for(--i, --j; j >=0; --i, --j)
  {
    mvec = frac * mt[j] + dobs * phi;
    C = frac * Ct[j].array();
    mvNormalDrawCov(mvec,C,1.0,p,phi);
    F = resp.block(0,j,n,p);
    yvec = resp.col(i);
    err = yvec - F*phi;
    Err[j] = err;
    vPhi[j] = phi;
  }
}

void stvar::mvNormalDrawCov(const VectorXd& mean, const MatrixXd& cor,
			    double sigma, int dim, VectorXd& out)
{
  LLTMd lltofcor(cor);
  if(lltofcor.info() != Eigen::Success)
    error("Cholesky failure for correlation matrix in multivariate normal sampling!");

  for(int i = 0; i != dim; ++i)
    out(i) = norm_rand();
  out.noalias() = lltofcor.matrixL()*out*sigma + mean;
}

void stvar::formSpx()
{
  int i, j;
  spx = MatrixXd::Zero(n,p);
  mvar = 0.0;
  for(i = 0; i != Tinf; ++i)
  {
    spx += Ft[i].cwiseProduct(vPhi[i].transpose().replicate(n,1));
    mvar += vVar[i];
  }
  spx /= (double)Tinf;
  mvar /= (double)Tinf;
  for(i = 0; i != p; ++i)
  {
    spresid[i] = spresp;
    for(j = 0; j != i; ++j)
      spresid[i] -= spx.col(j).cwiseProduct(rho.col(j));
    for(j = i+1; j != p; ++j)
      spresid[i] -= spx.col(j).cwiseProduct(rho.col(j));
  }
}
void stvar::drawRho()
{
  MatrixXd Fp, solFp, prec, umat, umatprec;
  VectorXd  mean, resid, coefi, newcoefi;
  VectorXd newrho(n);
  LLTMd lltofprec;
  int i, j;
  for(i = 0; i != p; ++i)
  {
    Fp = spx.col(i).asDiagonal();
    solFp = lltcorrRes.solve(Fp)/mvar;

    umat = lltRho[i].matrixU();
    ViewUmd uview = umat.triangularView<Upper>();
    ViewLcmd lview = lltRho[i].matrixL();
    
    prec = uview*Fp*solFp*lview;
    prec += MatrixXd::Identity(n,n);
    
    lltofprec.compute(prec);
    if(lltofprec.info() != Eigen::Success)
      error("Cholesky failure while drawing rho!");
    mean = solFp.transpose()*spresid[i];
    mean.noalias() = lview * lltofprec.solve(uview*mean);
    
    for(j = 0; j != n; ++j)
      newrho(j) = norm_rand();
    
    umatprec = lltofprec.matrixU();
    ViewUmd uviewprec = umatprec.triangularView<Upper>();
    newrho = lview * uviewprec.solve(newrho);
    newrho += mean;
    
    coefi = spx.col(i).cwiseProduct(rho.col(i));
    newcoefi = spx.col(i).cwiseProduct(newrho);
    for(j = 0; j != i; ++j)	// evade in loop if
    {
      resid = spresid[i] + coefi - newcoefi;
      spresid[i] = resid;
    }
    for(j = i+1; j != p; ++j);
    {
      resid = spresid[i] + coefi - newcoefi;
      spresid[i] = resid;
    }
    loglikRho[i] = evalLogLikRho(newrho,i);
    rho.col(i) = newrho;
  }
  for(i = 0; i != Tinf; ++i)
    Ft[i] = resp.block(0, i, n, p).cwiseProduct(rho);
}

void stvar::drawThetaRho()
{
  VectorXd newtheta, rhoi;
  double newthetaj, othetaj, newloglik;
  double logratio;
  MatrixXd newcorr(n,n);
  LLTMd lltofCorr;
  double u;
  for(int i = 0; i != p; ++i)
  {
    newtheta = thetaRho[i];
    rhoi = rho.col(i);
    for(int j = 0; j != d; ++j)
    {
      othetaj = newtheta(j);
      newthetaj = pKerThetaRho -> draw(othetaj);
      newtheta(j) = newthetaj;
      newloglik = stvarBase::evalLogLikRho(newtheta, rhoi, newcorr, lltofCorr);
      logratio = newloglik - loglikRho[i] + priThetaRho -> logpdf(newthetaj)
	- priThetaRho -> logpdf(othetaj);
      logratio += pKerThetaRho -> logdensity(newthetaj, othetaj)
	- pKerThetaRho -> logdensity(othetaj, newthetaj);
      u = runif(0.0, 1.0);
      if(log(u)<logratio)
      {
	thetaRho[i](j) = newthetaj;
	loglikRho[i] = newloglik;
	lltRho[i] = lltofCorr;
      }
    }
  }
}
double stvar::evalLogLikRho(const VectorXd& rho, int i)
{
  MatrixXd lmat;
  VectorXd solrho;
  double loglik, quadrho;
  LLTMd lltofCorr = lltRho[i];
  solrho = lltofCorr.solve(rho);
  quadrho = rho.dot(solrho);
  lmat = lltofCorr.matrixL();
  loglik = -lmat.diagonal().array().log().sum() - 0.5*quadrho;
  return loglik;
}
void stvar::mcmc()
{
  int diff = 1 - nburn;
  for(int i = 1; i <= nmcmc; ++i, ++diff)
  {
    forwardFilter();
    drawVar();
    forwardSmooth();
    backwardSample();
    drawThetaRes();
    formSpx();
    drawRho();
    drawThetaRho();
    if(diff >= 0 && diff%nthin == 0)
    {
      hthetaRes.push_back(thetaRes);
      hlltcorrRes.push_back(lltcorrRes);
      hvPhi.push_back(vPhi);
      hthetaRho.push_back(thetaRho);
      hlltRho.push_back(lltRho);
      hrho.push_back(rho);
    }
  }
}

void stvar::tvarmcmc()
{
  int diff = 1 - nburn;
  for(int i = 1; i <= nmcmc; ++i, ++diff)
  {
    tvarForwardFilter();
    drawVar();
    tvarForwardSmooth();
    tvarBackwardSample();
    drawThetaRes();
    if(diff >= 0 && diff%nthin == 0)
    {
      hthetaRes.push_back(thetaRes);
      hlltcorrRes.push_back(lltcorrRes);
      hvPhi.push_back(vPhi);
    }
  }
}

void stvar::predict(const MatrixXd& newx, const MatrixXd& newystart,
		    double* predout, double* newrhoout)
{
  int newn = newx.rows();
  MatrixXd newrho(newn,p);
  MatrixXd crosscorrrho(newn,n), crosscorrres(newn,n);
  MatrixXd pred(newn,T), newF, F;
  VectorXd err;
  double *curpred = predout, *currho = newrhoout;
  int bpred = newn*T, brho = newn*p;
  pred.leftCols(p) = newystart;
  for(int i = 0; i != nrec; ++i, curpred += bpred, currho += brho)
  {
    for(int j = 0; j != p; ++j)
    {
      pcorrRho -> evalCrossMat(newx,xmat,hthetaRho[i][j],crosscorrrho);
      newrho.col(j) = crosscorrrho * hlltRho[i][j].solve(hrho[i].col(j));
    }
    pcorrRes -> evalCrossMat(newx,xmat,hthetaRes[i],crosscorrres);
    for(int j = 0, k = p; j != Tinf; ++j, ++k)
    {
      F = resp.block(0,j,n,p);//.cwiseProduct(hrho[i]);
      newF = pred.block(0,j,newn,p);//.cwiseProduct(newrho);
      F = F.cwiseProduct(hrho[i]);
      newF = newF.cwiseProduct(newrho);
      pred.col(k) = newF*hvPhi[i][j];
      err = resp.col(k) - F * hvPhi[i][j];
      pred.col(k) += crosscorrres * hlltcorrRes[i].solve(err);
    }
    std::copy(pred.data(), pred.data()+bpred, curpred);
    std::copy(newrho.data(), newrho.data()+brho, currho);
  }
}

void stvar::tvarpredict(const MatrixXd& newx, const MatrixXd& newystart,
			double* predout)
{
  int newn = newx.rows();
  MatrixXd crosscorrres(newn,n);
  MatrixXd pred(newn,T), newF, F;
  VectorXd err;
  double* curpred = predout;
  int bpred = newn*T;
  pred.leftCols(p) = newystart;
  for(int i = 0; i != nrec; ++i, curpred += bpred)
  {
    pcorrRes -> evalCrossMat(newx,xmat,hthetaRes[i],crosscorrres);
    for(int j = 0, k = p; j != Tinf; ++j, ++k)
    {
      F = resp.block(0,j,n,p);
      newF = pred.block(0,j,newn,p);
      pred.col(k) = newF*hvPhi[i][j];
      err = resp.col(k) - F * hvPhi[i][j];
      pred.col(k) += crosscorrres * hlltcorrRes[i].solve(err);
    }
    std::copy(pred.data(), pred.data()+bpred, curpred);
  }
}
