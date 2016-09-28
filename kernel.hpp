#ifndef KERNEL_H
#define KERNEL_H
#include<cmath>
#include<cassert>
#include<R.h>
#include<Rmath.h>
#define MAX(a,b) ((a)>(b)? (a):(b))
#define MIN(a,b) ((a)<(b)? (a):(b))
class kernel
{
public:
    virtual double draw(double)const = 0;
    virtual double density(double,double)const=0;
    virtual double logdensity(double,double)const=0;
    virtual ~kernel(){};
};
// uniform random walk kernel assume to be symetric
class uniformKernel: public kernel
{
 public:
    uniformKernel(double step_in, double ub_in, double lb_in, bool isShift_in)
    {
      assert((isShift_in && step_in>0.0) || (!isShift_in && step_in>1.0));
      assert(ub_in>lb_in);
      step = step_in;
      isShift = isShift_in;
      ub = ub_in;
      lb = lb_in;
    }
    double draw(double center) const
    {
      double rup = isShift? (center+step):(center*step);
      double rlow = isShift? (center-step):(center/step);
	
      rup = MIN(rup,ub);
      rlow = MAX(rlow,lb);
      double res = runif(rlow,rup);
      return res;
    }
    double density(double from, double to) const
    {
      return 1.0;
    }
    double logdensity(double from, double to) const
    {
      return 0.0;
    }
private:
    double step;
    double ub,lb;
    bool isShift; 		// shift or shrink 
};

// normal random walk kernel
class normalKernel: public kernel
{
public:
    normalKernel(double sigma_in)
    {
	assert(sigma_in>0.0);
	sigma = sigma_in;
    }
    double draw(double center) const 
    {
	return rnorm(center,sigma);
    }
    double density(double from, double to) const 
    {
	return dnorm(to, from, sigma, 0);
    }
    double logdensity(double from, double to) const 
    {
	return dnorm(to, from, sigma, 1);
    }
private:
    double sigma;
};

class lognormalKernel: public kernel
{
public:
    lognormalKernel(double sigma_in)
    {
	assert(sigma_in>0.0);
	sigma = sigma_in;
    }
    double draw(double center) const
    {
	return rlnorm(log(center),sigma);
    }
    double density(double from, double to) const
    {
	return dlnorm(to, log(from), sigma, 0);
    }
    double logdensity(double from, double to) const
    {
	return dlnorm(to, log(from), sigma, 1);
    }
private:
    double sigma; 		// the standard deviation of log(X)
};
#endif
