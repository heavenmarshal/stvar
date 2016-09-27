#ifndef STVAR_H
#define STVAR_H
#include"stvarBase.hpp"
class stvar: public stvarBase
{
public:
  stvar(Mmat xmat_, Mmat resp_, int p_, const VectorXd& delta,
	double s0_, double n0_, int nmcmc_, int nburn_,
	int nthin_, int nrec_, double powrho, double powres,
	int priRho, const VectorXd& priParRho,
	int priRes, const VectorXd& priParRes,
	int kerRho, const VectorXd& parKerRho,
	int kerRes, const VectorXd& parKerRes,
	const VectorXd& thetaRes0, const VectorXd& thetaRho0,
	Mmat C0_, Mvec m0_):
    stvarBase(xmat_,resp_,p_,delta,s0_,n0_,nmcmc_,nburn_,
	      nthin_, nrec_, powrho, powres,priRho, priParRho,
	      priRes, priParRes, kerRho, parKerRho,
	      kerRes, parKerRes, thetaRes0),
    C0(C0_.data(),C0_.rows(),C0_.cols()),
    m0(m0_.data(),m0_.size()),
    loglikRho(p_,0.0)
  {
    rho = MatrixXd::Constant(n,p,1.0);
    spresid.resize(p);
    vPhi.resize(Tinf);
    mt.resize(Tinf);
    Ct.resize(Tinf);
    for(int i = 0; i != Tinf; ++i)
      Ft.push_back(resp.block(0, i, n, p));
    initThetaRho(thetaRho0);
  }
  
  stvar(Mmat xmat_, Mmat resp_, int p_, const VectorXd& delta,
	double s0_, double n0_, int nmcmc_, int nburn_,
	int nthin_, int nrec_, double powres,
	int priRes, const VectorXd& priParRes,
	int kerRes, const VectorXd& parKerRes,
	const VectorXd& thetaRes0,
	Mmat C0_, Mvec m0_):
    stvarBase(xmat_,resp_,p_,delta,s0_,n0_,nmcmc_,nburn_,
	      nthin_, nrec_, powres, priRes, priParRes,
	      kerRes, parKerRes, thetaRes0),
    C0(C0_.data(),C0_.rows(),C0_.cols()),
    m0(m0_.data(),m0_.size())
  {
    vPhi.resize(Tinf);
    mt.resize(Tinf);
    Ct.resize(Tinf);
  }

  void mcmc();
  void tvarmcmc();
  void predict(const MatrixXd&, const MatrixXd&, double*, double*);
  void tvarpredict(const MatrixXd&, const MatrixXd&, double*);
private:
  
  Mmat C0;
  Mvec m0;
  MatrixXd rho;
  vector<VectorXd> vPhi;
  MatrixXd spx;
  vector<VectorXd> thetaRho, spresid;
  vector<MatrixXd> precRho, hrho;
  vector<double> loglikRho;
  vvVecd hthetaRho, hvPhi;
  vvMatd hprecRho;
  vector<VectorXd> mt;
  vector<MatrixXd> Ct, Ft;
  void initThetaRho(const VectorXd&);
  void drawRho();
  void drawThetaRho();
  void formSpx();
  void forwardFilter();
  void tvarForwardFilter();
  void forwardSmooth();
  void tvarForwardSmooth();
  void backwardSample();
  void tvarBackwardSample();
  //double evalLogLikRho(const VectorXd&, int);
  static void mvNormalDrawCov(const VectorXd&, const MatrixXd&, double, int, VectorXd&);
};

#endif
