#ifndef CORRELATION_H
#define CORRELATION_H
#include<cmath>
#include<Rmath.h>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<cassert>
using Eigen::MatrixXd;
using Eigen::VectorXd;

class correlation
{
public:
  virtual void evalCorrMat(const MatrixXd&, const VectorXd&, MatrixXd&) const = 0;
  virtual void evalCrossMat(const MatrixXd&, const MatrixXd&, const VectorXd&, MatrixXd&) const=0;
  virtual ~correlation() {}
};

class powerexpCorr: public correlation
{
public:
  powerexpCorr(int dim_in, VectorXd power_in)
  {
    assert(dim_in>0);
    for(int i=1; i != dim_in; ++i)
      assert(power_in(i)>0 && power_in(i)<=2);
    dim = dim_in;
    power = power_in;
  }
  void evalCorrMat(const MatrixXd&, const VectorXd&, MatrixXd&) const;
  void evalCrossMat(const MatrixXd&, const MatrixXd&, const VectorXd&, MatrixXd&) const;
private:
  int dim;
  VectorXd power;
};

#endif
