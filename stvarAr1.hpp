#ifndef STVARAR1_H
#define STVARAR1_H
#include"stvarBase.hpp"
class stvarAr1: public stvarBase
{
public:
  stvarAr1(Mmat xmat_, Mmat resp_, const VectorXd& delta,
	   double s0_, double n0_, int nmcmc_, int nburn_,
	   int nthin_, int nrec_, double powrho, double powres,
	   int priRho, const VectorXd& priParRho,
	   int priRes, const VectorXd& priParRes,
	   int kerRho, const VectorXd& parKerRho,
	   int kerRes, const VectorXd& parKerRes,
	   const VectorXd& thetaRes0, const VectorXd& thetaRho0,
	   double C0_, double m0_):
    stvarBase(xmat_,resp_,1,delta,s0_,n0_,nmcmc_,nburn_,
	      nthin_, nrec_, powrho, powres, priRho,priParRho,
	      priRes,priParRes,kerRho,parKerRho,
	      kerRes, parKerRes, thetaRes0),
    C0(C0_), m0(m0_)
  {
    rho = VectorXd::Constant(n,1.0);
    vPhi.resize(Tinf);
    mt.resize(Tinf), Ct.resize(Tinf);
    for(int i = 0; i != Tinf; ++i)
      Ft.push_back(resp.col(i));
    initThetaRho(thetaRho0);
  }
  stvarAr1(Mmat xmat_, Mmat resp_, const VectorXd& delta,
	   double s0_, double n0_, int nmcmc_, int nburn_,
	   int nthin_, int nrec_, double powres,
	   int priRes, const VectorXd& priParRes,
	   int kerRes, const VectorXd& parKerRes,
	   const VectorXd& thetaRes0,
	   double C0_, double m0_):
    stvarBase(xmat_,resp_,1,delta,s0_,n0_,nmcmc_,nburn_,
	      nthin_, nrec_, powres, priRes, priParRes,
	      kerRes, parKerRes, thetaRes0),
    C0(C0_), m0(m0_)
  {
    vPhi.resize(Tinf);
    mt.resize(Tinf), Ct.resize(Tinf);
  }

  void mcmc();
  void tvarmcmc();
  void predict(const MatrixXd&, const MatrixXd&, double*, double*);
  void tvarpredict(const MatrixXd&, const MatrixXd&, double*);
  void getThetaRho(double*){};	// temporary
  void getPhi(double*){};	// temporary
private:
  double C0, m0;
  VectorXd rho;
  vector<double> vPhi;
  VectorXd spx;
  VectorXd thetaRho;
  vector<VectorXd> hthetaRho;

  MatrixXd precRho;
  double loglikRho;
  vector<VectorXd> hrho;
  vector<MatrixXd> hprecRho;
  vvd hvPhi;
  vector<double> mt, Ct;
  vector<VectorXd> Ft;
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
};
#endif
