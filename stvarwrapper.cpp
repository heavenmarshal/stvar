#include"stvar.hpp"
#include"stvarAr1.hpp"
#include<algorithm>
extern "C"
{
void stvar_R(const double* xmat_, const double* resp_, const int* n_, const int* T_,
	     const int *d_, const int *p_, const double *delta_, const double* s0_,
	     const double *n0_, const int* nmcmc_, const int *nburn_, const int *nthin_,
	     const int *nrec_, const double *powrho_, const double *powres_, const int *priRho_,
	     const double *priParRho_, const int* lenPriRho_, const int* priRes_,
	     const double* priParRes_, const int* lenPriRes_, const int *kerRho_,
	     const double *parKerRho_, const int* lenKerRho_, const int* kerRes_,
	     const double *parKerRes_, const int* lenKerRes_,
	     const double *thetaRes0_, const double *thetaRho0_,
	     const double* C0_, const double *m0_, const double *newx_,
	     const double* newystart_, const int* newn_, double* predout_, double* rho_)
  {
    int n = *n_, T = *T_, d = *d_, p = *p_;
    Mmat xmat(xmat_, n, d);
    Mmat resp(resp_, n, T);
    Mvec delta(delta_,2);
    Mvec priParRho(priParRho_,*lenPriRho_);
    Mvec priParRes(priParRes_,*lenPriRes_);
    Mvec parKerRho(parKerRho_,*lenKerRho_);
    Mvec parKerRes(parKerRes_,*lenKerRes_);
    Mvec thetaRes0(thetaRes0_,d);
    Mvec thetaRho0(thetaRho0_,d);
    stvarBase* pstvar;
    if(p == 1)
      pstvar = new stvarAr1(xmat,resp,delta,*s0_,*n0_,*nmcmc_,*nburn_,*nthin_,*nrec_,
			    *powrho_,*powres_, *priRho_, priParRho, *priRes_,
			    priParRes, *kerRho_, parKerRho, *kerRes_,parKerRes,thetaRes0,
			    thetaRho0,*C0_,*m0_);
    else
    {
      Mmat C0(C0_,p,p);
      Mvec m0(m0_,p);
      pstvar = new stvar(xmat,resp,p,delta,*s0_,*n0_,*nmcmc_,*nburn_,*nthin_,*nrec_,
			 *powrho_,*powres_, *priRho_, priParRho, *priRes_,
			 priParRes, *kerRho_, parKerRho, *kerRes_,parKerRes,
			 thetaRes0, thetaRho0, C0,m0);
    }
    pstvar -> mcmc();
    Mmat newx(newx_,*newn_,d);
    Mmat newystart(newystart_, *newn_, p);
    vector<MatrixXd> predout;
    pstvar -> predict(newx,newystart,predout_,rho_);
  }
  void tvar_R(const double* xmat_, const double* resp_, const int* n_, const int* T_,
	      const int *d_, const int *p_, const double *delta_, const double* s0_,
	      const double *n0_, const int* nmcmc_, const int *nburn_, const int *nthin_,
	      const int *nrec_, const double *powres_, const int* priRes_,
	      const double* priParRes_, const int* lenPriRes_, const int* kerRes_,
	      const double *parKerRes_, const int* lenKerRes_, const double *thetaRes0_,
	      const double* C0_, const double *m0_, const double *newx_,
	      const double* newystart_, const int* newn_, double* predout_)
  {
    int n = *n_, T = *T_, d = *d_, p = *p_;
    Mmat xmat(xmat_, n, d);
    Mmat resp(resp_, n, T);
    Mvec delta(delta_,2);
    Mvec priParRes(priParRes_,*lenPriRes_);
    Mvec parKerRes(parKerRes_,*lenKerRes_);
    Mvec thetaRes0(thetaRes0_,d);
    stvarBase* pstvar;
    if(p == 1)
      pstvar = new stvarAr1(xmat,resp,delta,*s0_,*n0_,*nmcmc_,*nburn_,*nthin_,*nrec_,
			    *powres_,  *priRes_, priParRes, *kerRes_,parKerRes,thetaRes0,
			    *C0_,*m0_);
    else
    {
      Mmat C0(C0_,p,p);
      Mvec m0(m0_,p);
      pstvar = new stvar(xmat,resp,p,delta,*s0_,*n0_,*nmcmc_,*nburn_,*nthin_,*nrec_,
			 *powres_, *priRes_, priParRes, *kerRes_, parKerRes,
			 thetaRes0, C0,m0);
    }
    pstvar -> tvarmcmc();
    Mmat newx(newx_,*newn_,d);
    Mmat newystart(newystart_, *newn_, p);
    vector<MatrixXd> predout;
    pstvar -> tvarpredict(newx,newystart,predout_);
  }
  
}
