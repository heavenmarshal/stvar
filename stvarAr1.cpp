#include"stvarAr1.hpp"
#include<cmath>
#include<R.h>
#include<Rmath.h>
#include<algorithm>
void stvarAr1::initThetaRho(const VectorXd& theta0)
{
  MatrixXd corrrho(n,n);
  pcorrRho -> evalCorrMat(xmat,theta0,corrrho);
  LLTMd lltcorrrho(corrrho);
  if(lltcorrrho.info() != Eigen::Success)
    error("Cholesky failure while initializing correlation for rho!");
  thetaRho = theta0;
  precRho = lltcorrrho.solve(MatrixXd::Identity(n,n));
}
void stvarAr1::forwardFilter()
{
  double R, Rinv, Rinfo;
  double finfo, fqinvf, C = C0, mvec = m0;
  double siter = s0, niter = n0, diter;
  double quade, quadQe, bilFe;
  VectorXd solF, F;
  VectorXd yvec, err, sole; 
  for(int i = 1, j=0; i != T; ++i,++j)
  {
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    solF = lltcorrRes.solve(F)/siter;
    sole = lltcorrRes.solve(err)/siter;
    quade = err.dot(sole);
    bilFe = solF.dot(err);
    finfo = F.dot(solF);
    Rinv = 1.0/R;
    Rinv = Rinv + finfo;
    quadQe = bilFe*bilFe/Rinv;
    Rinfo = finfo/Rinv;
    fqinvf = finfo - finfo*Rinfo;
    C = R - R*fqinvf*R;
    mvec = mvec + R * bilFe - R * Rinfo * bilFe;
    diter = dobs * niter * siter;
    diter += siter * (quade-quadQe);
    niter = dobs * niter + n;
    siter = diter/niter;
    Ct[j] = C;
    mt[j] = mvec;
    nt[j] = niter;
    st[j] = siter;
  }
}
void stvarAr1::tvarForwardFilter()
{
  double R, Rinv, Rinfo;
  double finfo, fqinvf, C = C0, mvec = m0;
  double siter = s0, niter = n0, diter;
  double quade, quadQe, bilFe;
  VectorXd solF, F;
  VectorXd yvec, err, sole; 
  for(int i = 1, j=0; i != T; ++i,++j)
  {
    F = resp.col(j);
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    solF = lltcorrRes.solve(F)/siter;
    sole = lltcorrRes.solve(err)/siter;
    quade = err.dot(sole);
    bilFe = solF.dot(err);
    finfo = F.dot(solF);
    Rinv = 1.0/R;
    Rinv = Rinv + finfo;
    quadQe = bilFe*bilFe/Rinv;
    Rinfo = finfo/Rinv;
    fqinvf = finfo - finfo*Rinfo;
    C = R - R*fqinvf*R;
    mvec = mvec + R * bilFe - R * Rinfo * bilFe;
    diter = dobs * niter * siter;
    diter += siter * (quade-quadQe);
    niter = dobs * niter + n;
    siter = diter/niter;
    Ct[j] = C;
    mt[j] = mvec;
    nt[j] = niter;
    st[j] = siter;
  }
}

void stvarAr1::forwardSmooth()
{
  double R, Rinv, Rinfo;
  double finfo, fqinvf, C = C0, mvec = m0;
  double bilFe;
  VectorXd F, solF, yvec, err;

  for(int i = p, j = 0; i != T; ++i, ++j)
  {
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    solF = lltcorrRes.solve(F)/vVar[j];
    bilFe = solF.dot(err);
    finfo = F.dot(solF);
    Rinv = 1.0/R;
    Rinv += finfo;
    Rinfo = finfo/Rinv;
    fqinvf = finfo - finfo*Rinfo;

    mvec += R*bilFe - R*Rinfo*bilFe;
    C = R - R*fqinvf*R;
    mt[j] = mvec;
    Ct[j] = C;
  }
}

void stvarAr1::tvarForwardSmooth()
{
  double R, Rinv, Rinfo;
  double finfo, fqinvf, C = C0, mvec = m0;
  double bilFe;
  VectorXd F, solF, yvec, err;

  for(int i = p, j = 0; i != T; ++i, ++j)
  {
    F = resp.col(j);
    yvec = resp.col(i);
    err = yvec - F*mvec;
    R = C/dsys;
    solF = lltcorrRes.solve(F)/vVar[j];
    bilFe = solF.dot(err);
    finfo = F.dot(solF);
    Rinv = 1.0/R;
    Rinv += finfo;
    Rinfo = finfo/Rinv;
    fqinvf = finfo - finfo*Rinfo;

    mvec += R*bilFe - R*Rinfo*bilFe;
    C = R - R*fqinvf*R;
    mt[j] = mvec;
    Ct[j] = C;
  }
}

void stvarAr1::backwardSample()
{
  int i = T-1, j = Tinf-1;
  double phi, mvec = mt[j], C = Ct[j];
  double frac = 1.0-dobs;
  VectorXd err, yvec = resp.col(i),  F = Ft[j];
  phi = rnorm(mvec,C);
  err = yvec - F*phi;
  Err[j] = err;
  vPhi[j] = phi;
  for(--i, --j; j >=0; --i, --j)
  {
    mvec = frac * mt[j] + dobs * phi;
    C = frac * Ct[j];
    phi = rnorm(mvec,C);
    F = Ft[j];
    yvec = resp.col(i);
    err = yvec - F*phi;
    Err[j] = err;
    vPhi[j] = phi;
  }
}

void stvarAr1::tvarBackwardSample()
{
  int i = T-1, j = Tinf-1;
  double phi, mvec = mt[j], C = Ct[j];
  double frac = 1.0-dobs;
  VectorXd err, yvec = resp.col(i), F = resp.col(j);
  phi = rnorm(mvec,C);
  err = yvec - F*phi;
  Err[j] = err;
  vPhi[j] = phi;
  for(--i, --j; j >=0; --i, --j)
  {
    mvec = frac * mt[j] + dobs * phi;
    C = frac * Ct[j];
    phi = rnorm(mvec,C);
    F = resp.col(j);
    yvec = resp.col(i);
    err = yvec - F*phi;
    Err[j] = err;
    vPhi[j] = phi;
  }
}

void stvarAr1::formSpx()
{
  spx = VectorXd::Zero(n);
  mvar = 0.0;
  for(int i = 0; i != Tinf; ++i)
  {
    spx += Ft[i]*vPhi[i];
    mvar += vVar[i];
  }
  spx /= (double)Tinf;
  mvar /= (double)Tinf;
}
void stvarAr1::drawRho()
{
  MatrixXd Fp, solFp, prec, corrrho(n,n);
  VectorXd likmean;
  LLTMd lltofrho;
  Fp = spx.asDiagonal();
  solFp = lltcorrRes.solve(Fp)/mvar;
  prec = Fp*solFp+precRho;
  likmean = solFp.transpose()*spresp;
  mvNormalDrawPrec(likmean,prec,1.0,n,rho);
  loglikRho = evalLogLikRho(thetaRho,rho,corrrho,lltofrho);
  precRho = lltofrho.solve(MatrixXd::Identity(n,n));
}

void stvarAr1::drawThetaRho()
{
  VectorXd newtheta;
  double newthetaj, othetaj, newloglik;
  double logratio, u;
  MatrixXd newcorr(n,n);
  LLTMd lltofCorr;
  newtheta = thetaRho;
  for(int j = 0; j != d; ++j)
  {
    othetaj = newtheta(j);
    newthetaj = pKerThetaRho -> draw(newtheta(j));
    newtheta(j) = newthetaj;
    newloglik = evalLogLikRho(newtheta, rho, newcorr, lltofCorr);
    logratio = newloglik - loglikRho + priThetaRho -> logpdf(newthetaj)
	- priThetaRho -> logpdf(othetaj);
    logratio += pKerThetaRho -> logdensity(newthetaj, othetaj)
      - pKerThetaRho -> logdensity(othetaj, newthetaj);
    u = runif(0.0, 1.0);
    if(log(u)<logratio)
    {
      thetaRho(j) = newthetaj;
      loglikRho = newloglik;
      precRho = lltofCorr.solve(MatrixXd::Identity(n,n));
    }
  }
}
void stvarAr1::mcmc()
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
      hthetaRho.push_back(thetaRho);
      hprecRho.push_back(precRho);
      hlltcorrRes.push_back(lltcorrRes);
      hrho.push_back(rho);
      hvPhi.push_back(vPhi);
    }
  }
}
void stvarAr1::tvarmcmc()
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
void stvarAr1::predict(const MatrixXd& newx, const MatrixXd& newystart,
		       double* predout, double* newrhoout)
{
  int newn = newx.rows();
  VectorXd newrho(newn);
  MatrixXd crosscorrrho(newn,n), crosscorrres(newn,n);
  MatrixXd pred(newn,T), newF, F;
  VectorXd err;

  double *curpred = predout, *currho = newrhoout;
  int bpred = newn*T, brho = newn;
  pred.leftCols(p) = newystart;
  for(int i = 0; i != nrec; ++i, curpred+=bpred, currho += brho)
  {
    pcorrRho -> evalCrossMat(newx,xmat,hthetaRho[i],crosscorrrho);
    newrho = crosscorrrho * hprecRho[i] * hrho[i];
    pcorrRes -> evalCrossMat(newx,xmat,hthetaRes[i],crosscorrres);
    for(int j = 0, k = p; j != Tinf; ++j, ++k)
    {
      F = resp.col(j).cwiseProduct(hrho[i]);
      err = resp.col(k) - F * hvPhi[i][j];
      newF = pred.col(j).cwiseProduct(newrho);
      pred.col(k) = newF*hvPhi[i][j];
      pred.col(k) += crosscorrres * hlltcorrRes[i].solve(err);
    }
    std::copy(pred.data(), pred.data()+bpred, curpred);
    std::copy(newrho.data(), newrho.data()+brho, currho);
  }
}
void stvarAr1::tvarpredict(const MatrixXd& newx, const MatrixXd& newystart,
			   double* predout)
{
  int newn = newx.rows();
  MatrixXd crosscorrres(newn,n);
  MatrixXd pred(newn,T), newF, F;
  VectorXd err;
  double *curpred = predout;
  int bpred = newn*T;
  pred.leftCols(p) = newystart;
  for(int i = 0; i != nrec; ++i, curpred += bpred)
  {
    for(int j = 0, k = p; j != Tinf; ++j, ++k)
    {
      F = resp.col(j);
      err = resp.col(k) - F * hvPhi[i][j];
      newF = pred.col(j);
      pred.col(k) = newF*hvPhi[i][j];
      pred.col(k) += crosscorrres * hlltcorrRes[i].solve(err);
    }
    std::copy(pred.data(), pred.data()+bpred, curpred);
  }
}
