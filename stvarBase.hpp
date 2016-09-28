#ifndef STVARBASE_H
#define STVARBASE_H
#include<Eigen/Dense>
#include<Eigen/Core>
#include<vector>
#include"kernel.hpp"
#include"correlation.h"
#include"distb.h"

#define GAMMA   101
#define INVGA   102

#define UNIF    201
#define LNORMAL 202
#define NORMAL  203

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::LLT;
using std::vector;
typedef Map<const MatrixXd> Mmat;
typedef Map<const VectorXd> Mvec;
typedef LLT<MatrixXd> LLTMd;
typedef vector<vector<double> > vvd;
typedef vector<vector<VectorXd> > vvVecd;
typedef vector<vector<MatrixXd> > vvMatd;
typedef vector<vector<LLTMd> > vvLLTMd;

class stvarBase
{
public:
  stvarBase(Mmat, Mmat, int, const VectorXd&, double, double, int, int,
	    int, int, double, double, int, const VectorXd&, int, const VectorXd&,
	    int, const VectorXd&, int, const VectorXd&, const VectorXd&);
  stvarBase(Mmat, Mmat, int, const VectorXd&, double, double, int, int,
	    int, int, double, int, const VectorXd&, int, const VectorXd&,
	    const VectorXd&);
  virtual void mcmc(){};
  virtual void tvarmcmc(){};
  virtual void predict(const MatrixXd&, const MatrixXd&, double*, double*){};
  virtual void tvarpredict(const MatrixXd&, const MatrixXd&, double*){};
  virtual void getThetaRho(double*){};
  virtual void getPhi(double*){};
  virtual ~stvarBase()
    {
      delete pcorrRho;
      delete pcorrRes;
      delete priThetaRho;
      delete priThetaRes;
      delete pKerThetaRho;
      delete pKerThetaRes;
    }
protected:
  Mmat xmat, resp;
  int n, T, d, p, Tinf;
  double dobs, dsys;
  VectorXd spresp;
  
  VectorXd thetaRes;
  MatrixXd corrRes;
  LLTMd lltcorrRes;
  vector<double> vVar;
  double mvar, s0, n0;
  vector<double> st, nt;
  vector<VectorXd> Err;
  
  // historical data
  vector<VectorXd> hthetaRes;
  vector<LLTMd> hlltcorrRes;
  vector<MatrixXd> hcorrRes;
  // variables for mcmc sampling
  int nmcmc, nburn, nthin, nrec;
  correlation *pcorrRho, *pcorrRes;
  distb *priThetaRho, *priThetaRes;
  kernel *pKerThetaRho, *pKerThetaRes;

  double loglikRes, logPriRes;
  
  void setPriorRho(int, const VectorXd&);
  void setPriorRes(int, const VectorXd&);
  void setKernelRho(int, const VectorXd&);
  void setKernelRes(int, const VectorXd&);
  void initThetaRes(const VectorXd&);
  void evalLoglikRes();
  double evalNewLoglikRes(const VectorXd&, MatrixXd&, LLTMd&);
  void drawThetaRes();
  void drawVar();
  double evalLogLikRho(const VectorXd&, const VectorXd&, MatrixXd&, LLTMd&);
  static void mvNormalDrawPrec(const VectorXd&, const MatrixXd&, double, int, VectorXd&);
};

#endif
