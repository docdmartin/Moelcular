#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "Matrix3x3.h"
#include "Common.h"
#include "Kernel.h"

using namespace std;

class Solver {

public:

  Solver( size_t n, vector<Matrix3x3>& hessian, vector<double> frequency, double* electric_potential, Kernel& kernel );
  ~Solver();

  vector<pair<double, double>> CalculateResponse( double weight );

  void Reset();

  void Print();

private:

  bool  Solve             ();
  void  RealMatrixMultiply( double* input, double* output );
  void  iteration         ();

  vector<Matrix3x3>& mHessian;
  Kernel& mKernel;

  size_t  mLength;
  vector<double> mFrequencyVec;
  double  mFrequency, mFreqSq;
  double* mElectricPotential;

  double* intermediate_vec;
  double* real_response;
  double* imag_response;
  double* mAb;
  double* mb;
  double* solve_response;

  double threshold, mConvergenceThreshold;
  double alpha, beta, omega, rho, tt, ts;
  double *a, *b, *c, *d, *g, *f;
  double *r, *r_hat, *p, *v, *h, *s, *t;
  double *r_hat_end, *v_end, *t_end;

};

#endif
