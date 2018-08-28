#include <numeric>      // std::inner_product

#include "math/LinearAlgebra.h"

LinearAlgebra::LinearAlgebra(){
  mNumberNodes             = 0;
  mRealDimension           = 0;
  mSize                    = 0;
  mSolverConvergeThreshold = 0.0;
  mToleranceSq             = 0.001 * 0.001;
}

LinearAlgebra::~LinearAlgebra(){

}


void LinearAlgebra::SetDimension(int num_of_nodes) {
  mNumberNodes             = num_of_nodes;
  mRealDimension           = 3 * mNumberNodes;
  mSize                    = 2 * mRealDimension;

  B.reserve( mSize );
  r.reserve( mSize );
  z.reserve( mSize );
  p.reserve( mSize );
  w.reserve( mSize );

  mConstantTerm.reserve( mRealDimension );
  mSolution.reserve( mSize );
}


void LinearAlgebra::SetConstantVector(vector<double>& b) {
  copy(b.begin(),b.end(),back_inserter(mConstantTerm));
}

void LinearAlgebra::SetConnection(int node_1, int node_2, vector<double> dp, double k) {
  mHessianMatrix.push_back( {3*node_1, 3*node_2, {dp[0], dp[1], dp[2]}, k} );
}

void LinearAlgebra::GetPotential(double* real_x, double* imag_x) {
  *real_x = 0.0;
  *imag_x = 0.0;

  int imag_offset = mRealDimension;
  for (int index = 0; index < mRealDimension; ++index) {
    *real_x += mConstantTerm[index] * mSolution[index];
    *imag_x += mConstantTerm[index] * mSolution[index + imag_offset];
  }
}


void LinearAlgebra::PreProcess() {

  int imag_offset = mRealDimension;
  fill(B.begin(), B.begin()+imag_offset, 0.0);
  fill(mSolution.begin(), mSolution.begin(), 0.0);

  double sq_term;
  for(OuterProduct op : mHessianMatrix) {
      for(int cnt = 0; cnt < 3; ++cnt) {
        sq_term = op.mDeltaPosition[cnt] * op.mDeltaPosition[cnt];
        B[op.mNodeIndex1 + cnt] += sq_term;
        B[op.mNodeIndex2 + cnt] += sq_term;
      }
  }

  mSolverConvergeThreshold = 0.0;
  int index       = 0;
  // Find the inverse of A's true diagonal to improve speed
  for (int node = 0; node < mNumberNodes; ++node) {
    for (int diag = 0; diag < 3; ++diag) {
      // Set inverse, multiplied by 2 ... once for real and again for imaginary
      B[index              ] = 1.0 / B[index];
      B[index + imag_offset] = B[index];

      // Initialize solution to be constant vector
      // Leave imaginary part zero
      mSolution[index] = mConstantTerm[index];

      // imaginary b is always zero
      double C = B[index] * mConstantTerm[index];
      // Track if there is a delta x to solve for

      mSolverConvergeThreshold += C * C;
      ++index;
    }
  }

  mSolverConvergeThreshold *= mToleranceSq;
}


void LinearAlgebra::PerformCG() {

  // If b-vector is zero, then terminateSolverThreshold will be zero
  if (mSolverConvergeThreshold == 0.0) {
    fill(mSolution.begin(), mSolution.end(), 0.0f);
    // No solver iterations required, solution is trivial
    return;
  }

  size_t index;

  // Allows for initial guess of x to be non-zero
  // z = A * x
  multiplyMatrixA(mSolution, z);
  for (index = 0; index < mSize; ++index) {
    r[index] = mConstantTerm[index] - z[index];
    z[index] = B[index] * r[index];
    p[index] = z[index];
  }

  // beta = sum(r * z)
  double beta = inner_product(r.begin(), r.end(), z.begin(), 0.0f);

  for (size_t curr_iter = 0; curr_iter < mSize; ++curr_iter) {
    // w = A * p
    multiplyMatrixA(p, w);

    // alpha = beta / sum(w' * p)
    double alpha = inner_product(w.begin(), w.end(), p.begin(), 0.0f);
    alpha = beta / alpha;

    /**
     * x = x + alpha * p
     * r = r - alpha * w
     * z = B * r
     * beta = z' * r
     * gamma = beta / beta_old
     */
    double gamma = 1.0 / beta;
    beta = 0.0;
    double threshold = 0.0;
    for (index = 0; index < mSize; ++index) {
      mSolution[index] += alpha    * p[index];
      r        [index] -= alpha    * w[index];
      z        [index]  = B[index] * r[index];
      beta             += z[index] * r[index];
      threshold        += z[index] * z[index];
    }

    // Check to see if this has completed
    if (threshold < mSolverConvergeThreshold) break;

    gamma *= beta;
    // p = z + gamma * p
    for (index = 0; index < mSize; ++index)
      p[index] = z[index] + gamma * p[index];
  }
}

void LinearAlgebra::multiplyMatrixA(vector<double>& x, vector<double>& b) {
  /*
      [ H            -omega * I ] [X_real]  = [ b ]
      [ omega * I     H         ] [X_imag]
  */
//      double               mOmega;
//      vector<OuterProduct> mHessianMatrix;

    int imag_offset = mRealDimension;
    for(int cnt = 0; cnt < imag_offset; ++cnt) {
      b[cnt              ] = -mOmega * x[cnt + imag_offset];
      b[cnt + imag_offset] =  mOmega * x[cnt              ];
    }

    double dot_prod;
    double tmp;
    int    imag_index1;
    int    imag_index2;
    for(OuterProduct op : mHessianMatrix) {
        // real portion of X
        dot_prod = (x[op.mNodeIndex1  ] - x[op.mNodeIndex2  ]) * op.mDeltaPosition[0]
                 + (x[op.mNodeIndex1+1] - x[op.mNodeIndex2+1]) * op.mDeltaPosition[1]
                 + (x[op.mNodeIndex1+2] - x[op.mNodeIndex2+2]) * op.mDeltaPosition[2];
        dot_prod *= op.mScalar;

        for(int cnt = 0; cnt < 3; ++cnt) {
          tmp = op.mDeltaPosition[cnt] * dot_prod;

          b[op.mNodeIndex1+cnt] += tmp;
          b[op.mNodeIndex2+cnt] -= tmp;
        }

        // imag portion of X
        imag_index1 = op.mNodeIndex1 + imag_offset;
        imag_index2 = op.mNodeIndex2 + imag_offset;

        dot_prod = (x[imag_index1  ] - x[imag_index2  ]) * op.mDeltaPosition[0]
                 + (x[imag_index1+1] - x[imag_index2+1]) * op.mDeltaPosition[1]
                 + (x[imag_index1+2] - x[imag_index2+2]) * op.mDeltaPosition[2];
        dot_prod *= op.mScalar;

        for(int cnt = 0; cnt < 3; ++cnt) {
          tmp = op.mDeltaPosition[cnt] * dot_prod;

          b[imag_index1+cnt] += tmp;
          b[imag_index2+cnt] -= tmp;
        }
    }
}
