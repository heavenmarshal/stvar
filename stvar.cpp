#include"stvar.hpp"
#include<cmath>
#include<R.h>
#include<Rmath.h>
#include<algorithm>
#include<stdexcept>
#include<iostream>
void stvar::initThetaRho(const VectorXd& theta0)
{
  MatrixXd corrrho(n,n), preccorr;
  pcorrRho -> evalCorrMat(xmat,theta0,corrrho);
  LLTMd lltcorrrho(corrrho);
  if(lltcorrrho.info() != Eigen::Success)
    error("Cholesky failure while initializing correlation for rho!");
  preccorr = lltcorrrho.solve(MatrixXd::Identity(n,n));
  for(int i = 0; i != p; ++i)
  {
    thetaRho.push_back(theta0);
    precRho.push_back(preccorr);
  }
}

void stvar::forwardFilter()
{
  MatrixXd F, R, Rinv, Rinfo, solF;
  MatrixXd finfo, fqinvf, C = C0;
  VectorXd yvec, err, sole, bilFe, mvec = m0;
  double siter = s0, niter = n0, diter;
  double quade, quadQe;
  LLTMd lltrinv;
  for(int i = p, j=0; i != T; ++i,++j)
  {
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C / dsys;
    solF = lltcorrRes.solve(F)/siter;
    sole = lltcorrRes.solve(err)/siter;
    quade = err.dot(sole);
    bilFe = solF.transpose()*err;
    finfo = F.transpose()*solF;
    Rinv = R.ldlt().solve(MatrixXd::Identity(p,p));
    Rinv.noalias() = Rinv + finfo;
    lltrinv.compute(Rinv);
    if(lltrinv.info() == Eigen::NumericalIssue)
      error("Cholesky failure while filtering forward!");

    quadQe = bilFe.dot(lltrinv.solve(bilFe));
    Rinfo = lltrinv.solve(finfo);
    fqinvf = finfo - finfo*Rinfo;
    C = R - R*fqinvf*R;
    C /= siter;
    mvec = mvec + R*bilFe - R*Rinfo.transpose()*bilFe;
    diter = dobs * niter * siter;
    diter += siter*(quade-quadQe);
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
  MatrixXd F, Fl, R, Rinfo, solFl, lmat;
  MatrixXd finfo, fqinvf, C = C0;
  VectorXd yvec, err, bilFle, madd, mvec = m0;
  LLTMd lltrinfo;
  for(int i = p, j = 0; i != T; ++i, ++j)
  {
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    lmat = R.llt().matrixL();
    Fl = F * lmat;
    solFl = lltcorrRes.solve(Fl)/vVar[j];
    finfo = Fl.transpose() * solFl;
    bilFle = solFle.transpose() * err;
    Rinfo = finfo + MatrixXd::Identity(p,p);
    lltrinfo.compute(Rinfo);
    if(lltrinfo.info() != Eigen::NumericalIssue)
      error("Cholesky failure while smoothing forward!");

    madd = bilFle - finfo * lltrinfo.solve(bilFle);
    madd.noalias() = lmat * madd;
    mvec += madd;
    fqinvf = finfo * lltrinfo.solve(finfo);
    C = R - lmat * fqinvf * lmat.transpose();
    C.noalias() = 0.5*(C+C.transpose());
    mt[j] = mvec;
    Ct[j] = C;
  }
}
void stvar::tvarForwardFilter()
{
  MatrixXd F, R, Rinv, Rinfo, solF;
  MatrixXd finfo, fqinvf, C = C0;
  VectorXd yvec, err, sole, bilFe, mvec = m0;
  double siter = s0, niter = n0, diter;
  double quade, quadQe;
  LLTMd lltrinv;
  for(int i = p, j=0; i != T; ++i,++j)
  {
    F = resp.block(0,j,n,p);
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    solF = lltcorrRes.solve(F)/siter;
    sole = lltcorrRes.solve(err)/siter;
    quade = err.dot(sole);
    bilFe = solF.transpose()*err;
    finfo = F.transpose()*solF;
    Rinv = R.ldlt().solve(MatrixXd::Identity(p,p));
    Rinv.noalias() = Rinv + finfo;
    lltrinv.compute(Rinv);
    if(lltrinv.info() != Eigen::Success)
      error("Cholesky failure while filtering forward!");

    quadQe = bilFe.dot(lltrinv.solve(bilFe));
    Rinfo = lltrinv.solve(finfo);
    fqinvf = finfo - finfo*Rinfo;
    C = R - R*fqinvf*R;
    C /= siter;
    mvec = mvec + R*bilFe - R*Rinfo.transpose()*bilFe;
    diter = dobs * niter * siter;
    diter += siter*(quade-quadQe);
    niter = dobs * niter + n;
    siter = diter/niter;
    C.noalias() = 0.5*(C+C.transpose());
    C *= siter;
    Ct[j] = C;
    mt[j] = mvec;
    nt[j] = niter;
    st[j] = siter;
  }
}
void stvar::tvarForwardSmooth()
{
  MatrixXd F, R, Rinv, Rinfo, solF;
  MatrixXd finfo, fqinvf, C = C0;
  VectorXd yvec, err, bilFe, mvec = m0;
  LLTMd lltrinv;
  for(int i = p, j = 0; i != T; ++i, ++j)
  {
    F = resp.block(0,j,n,p);
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    solF = lltcorrRes.solve(F)/vVar[j];
    bilFe = solF.transpose()*err;
    finfo = F.transpose()*solF;
    Rinv = R.ldlt().solve(MatrixXd::Identity(p,p));
    Rinv.noalias() = Rinv + finfo;
    lltrinv.compute(Rinv);
    if(lltrinv.info() != Eigen::Success)
      error("Cholesky failure while smoothing forward!");
    Rinfo = lltrinv.solve(finfo);
    fqinvf = finfo - finfo*Rinfo;

    mvec += R*bilFe - R*Rinfo.transpose()*bilFe;
    C = R - R*fqinvf*R;
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
  try
  {
    mvNormalDrawCov(mvec, C, 1.0, p, phi);
  }
  catch(const std::runtime_error& e)
  {
    phi = mvec;
    warning("Cholesky failure while sampling backward, set to mean!");
  }
  err = yvec - F*phi;
  Err[j] = err;
  vPhi[j] = phi;
  for(--i, --j; j >=0; --i, --j)
  {
    mvec = frac * mt[j] + dobs * phi;
    C = frac * Ct[j];
    try
    {
      mvNormalDrawCov(mvec,C,1.0,p,phi);
    }
    catch(const std::runtime_error& e)
    {
      phi = mvec;
      warning("Cholesky failure while sampling backward, set to mean!");
    }
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
    C = frac * Ct[j];
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
    throw std::runtime_error("Cholesky error");
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
  MatrixXd Fp, prec, solFp, corrrho(n,n);
  VectorXd likmean, resid, coefi, newcoefi;
  VectorXd newrho(n);
  LLTMd lltofrho;
  int i, j;
  for(i = 0; i != p; ++i)
  {
    Fp = spx.col(i).asDiagonal();
    solFp = lltcorrRes.solve(Fp)/mvar;

    prec = Fp*solFp+precRho[i];
    likmean = solFp.transpose()*spresid[i];
    mvNormalDrawPrec(likmean,prec,1.0,n,newrho);

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
    loglikRho[i] = stvarBase::evalLogLikRho(thetaRho[i],newrho,corrrho,lltofrho);
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
	precRho[i] = lltofCorr.solve(MatrixXd::Identity(n,n));
      }
    }
  }
}
// double stvar::evalLogLikRho(const VectorXd& rho, int i)
// {
//   MatrixXd lmat;
//   VectorXd solrho;
//   double loglik, quadrho;
//   LLTMd lltofCorr = lltRho[i];
//   solrho = lltofCorr.solve(rho);
//   quadrho = rho.dot(solrho);
//   lmat = lltofCorr.matrixL();
//   loglik = -lmat.diagonal().array().log().sum() - 0.5*quadrho;
//   return loglik;
// }
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
      hprecRho.push_back(precRho);
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
      newrho.col(j) = crosscorrrho * hprecRho[i][j] * hrho[i].col(j);
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
