#include <cmath>

#include "math/LinearAlgebra.h"

LinearAlgebra::LinearAlgebra(){
  mNumberNodes             = 0;
  mRealDimension           = 0;
  mSize                    = 0;
  mSolverConvergeThreshold = 0.0;
  mToleranceSq             = 0.001 * 0.001;

  mNormalizedTranslation  = 0.0;
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

  mNormalizedTranslation = 1.0 / sqrt( static_cast<double>(mNumberNodes) );

  mRotationXY.reserve( mRealDimension );
  mRotationXZ.reserve( mRealDimension );
  mRotationYZ.reserve( mRealDimension );
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
  for(int cnt = 0; cnt < static_cast<int>(mSize); ++cnt) {
      if(cnt < imag_offset)
          B.push_back(-mOmega * mConstantTerm[cnt]);
      mSolution.push_back(0.0);
  }

  mSolverConvergeThreshold = 0.0;
  int index       = 0;
  for (int node = 0; node < mNumberNodes; ++node) {
    for (int diag = 0; diag < 3; ++diag) {

      mSolution[index] = B[index];

      mSolverConvergeThreshold += B[index] * B[index];
      ++index;
    }
  }

  mSolverConvergeThreshold *= mToleranceSq;
}


void LinearAlgebra::PerformCG() {

//    printMatrixA();
//    printMatrixB();

    int index;

    multiplyLinearSystem(mSolution, p);

    double gamma = 0.0;
    for (index = 0; index < mRealDimension; ++index) {
        r[index] = B[index] - z[index];
        p[index] = r[index];
        gamma   += r[index] * r[index];
    }

    double beta;
    int curr_iter;
    for ( curr_iter = 0; curr_iter < mRealDimension; ++curr_iter) {

        multiplyLinearSystem(p, w);

        double alpha = 0.0;
        for (index = 0; index < mRealDimension; ++index)
            alpha += w[index] * p[index];
        alpha = gamma / alpha;

        beta  = 1.0 / gamma;
        gamma = 0.0;

        for (index = 0; index < mRealDimension; ++index) {
            mSolution[index] += alpha    * p[index];
            r        [index] -= alpha    * w[index];
            gamma            += r[index] * r[index];
        }

        if (gamma < mSolverConvergeThreshold)
            break;

        cout << gamma << " > " << mSolverConvergeThreshold << endl;

        beta *= gamma;
        for (index = 0; index < mRealDimension; ++index)
            p[index] = r[index] + beta * p[index];
    }

    cout << gamma << " <= " << mSolverConvergeThreshold << endl;
    cout << curr_iter << " out of " << mRealDimension << " iterations were required" << endl;

    multiplyMatrixA(mSolution, p);
    beta = -1.0 / mOmega;
    for (index = 0; index < mRealDimension; ++index) {
        mSolution[index+mRealDimension] = mSolution[index];
        mSolution[index               ] = beta * p[index];
    }
}

void LinearAlgebra::multiplyLinearSystem(vector<double>& x, vector<double>& b) {
    double omega2 = mOmega * mOmega;
    multiplyMatrixA(x, z);
    multiplyMatrixA(z, b);
    for (int index = 0; index < mRealDimension; ++index)
        b[index] += omega2 * x[index];
}

void LinearAlgebra::multiplyMatrixA(vector<double>& x, vector<double>& b) {

    for(int cnt = 0; cnt < mRealDimension; ++cnt) {
      b[cnt] = 0.0;
    }

    double dot_prod;
    double tmp;
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
    }
}



void LinearAlgebra::printMatrixA(){

  cout << "Omega = " << mOmega << endl;
  cout << "A = zeros(" << mRealDimension << "," << mRealDimension << ");" << endl;

  for(OuterProduct op : mHessianMatrix) {

  cout << "B = " << op.mScalar << " * ["
       << op.mDeltaPosition[0] * op.mDeltaPosition[0] << ", "
       << op.mDeltaPosition[0] * op.mDeltaPosition[1] << ", "
       << op.mDeltaPosition[0] * op.mDeltaPosition[2] << ";"
       << op.mDeltaPosition[1] * op.mDeltaPosition[0] << ", "
       << op.mDeltaPosition[1] * op.mDeltaPosition[1] << ", "
       << op.mDeltaPosition[1] * op.mDeltaPosition[2] << ";"
       << op.mDeltaPosition[2] * op.mDeltaPosition[0] << ", "
       << op.mDeltaPosition[2] * op.mDeltaPosition[1] << ", "
       << op.mDeltaPosition[2] * op.mDeltaPosition[2] << "];" << endl;

  cout << "A("     << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
       << ") = A(" << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
       << ") + B;" << endl;
  cout << "A("     << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
       << ") = A(" << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
       << ") - B;" << endl;
  cout << "A("     << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
       << ") = A(" << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
       << ") - B;" << endl;
  cout << "A("     << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
       << ") = A(" << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
       << ") + B;" << endl;

     }
}

void LinearAlgebra::printMatrixB() {
  cout << "bb = [" << mConstantTerm[0];
  for (int index = 1; index < mRealDimension; ++index) {
      cout << ";" << mConstantTerm[index];
  }
  cout << "];" << endl;

  cout << "b = [" << B[0];
  for (int index = 1; index < mRealDimension; ++index) {
      cout << ";" << B[index];
  }
  cout << "];" << endl;

  cout << "C = A*A + Omega * Omega * eye(" << mRealDimension << ");" << endl;
  cout << "y = C\\b;" << endl;
  cout << "x = -1/Omega * A * y;" << endl;
}


















/*

void LinearAlgebra::PreProcess() {

for(double d : mConstantTerm)
cout << d << endl;
//  removeProjIntoZeroEigen( mConstantTerm );
for(double d : mConstantTerm)
cout << d << endl;

  int imag_offset = mRealDimension;
  for(int cnt = 0; cnt < static_cast<int>(mSize); ++cnt) {
      if(cnt < imag_offset)
          B.push_back(0.0);
      mSolution.push_back(0.0);
  }
//  fill(B.begin(), B.begin()+imag_offset, 0.0);
//  fill(mSolution.begin(), mSolution.begin()+mSize, 0.0);

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

      B[index              ] = 1.0;// / B[index];
      B[index + imag_offset] = B[index];

  cout << "B: " << B[index] << ", " << B[index + imag_offset] << endl;

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
*/

/*
void LinearAlgebra::PerformCG() {

  // If b-vector is zero, then terminateSolverThreshold will be zero
  if (mSolverConvergeThreshold == 0.0) {
    fill(mSolution.begin(), mSolution.end(), 0.0f);
    // No solver iterations required, solution is trivial
    return;
  }

cout << "mSolution size = " << mSolution.size() << endl;
  for(double d : mSolution)
  cout << "mSolution : " << d << endl;

  size_t index;

  // Allows for initial guess of x to be non-zero
  // z = A * x
  cout << "Calling first matrix multiply" << endl;
  multiplyMatrixA(mSolution, z);
  cout << "Completed first matrix multiply" << endl;

  double beta = 0.0; // beta = sum(r * z)
  for (index = 0; index < mSize; ++index) {
    if(index < static_cast<size_t>(mRealDimension))
      r[index] = mConstantTerm[index] - z[index];
    else
      r[index] = -z[index];
    z[index] = B[index] * r[index];
    p[index] = z[index];

    beta += r[index] * z[index];
  }

//  removeProjIntoZeroEigen( r );

  size_t curr_iter;
  for ( curr_iter = 0; curr_iter < mSize; ++curr_iter) {

    // w = A * p
    multiplyMatrixA(p, w);

    // alpha = beta / sum(w' * p)
    double alpha = 0.0;
    for (index = 0; index < mSize; ++index)
        alpha += w[index] * p[index];
    alpha = beta / alpha;


     // x = x + alpha * p
     // r = r - alpha * w
     // z = B * r
     // beta = z' * r
     // gamma = beta / beta_old

    double gamma = 1.0 / beta;
    beta = 0.0;
    double threshold = 0.0;
    for (index = 0; index < mSize; ++index) {
      mSolution[index] += alpha    * p[index];
      r        [index] -= alpha    * w[index];
      z        [index]  = B[index] * r[index]; // This is the issue B * r will map into kernal
      beta             += z[index] * r[index];
      threshold        += z[index] * z[index];
    }

    // Check to see if this has completed
    if (threshold < mSolverConvergeThreshold)
        break;

//    removeProjIntoZeroEigen( r );

    cout << threshold << " > " << mSolverConvergeThreshold << endl;

    gamma *= beta;
    // p = z + gamma * p
    for (index = 0; index < mSize; ++index)
      p[index] = z[index] + gamma * p[index];
  }

  cout << curr_iter << " out of " << mSize << " iterations were required" << endl;
}
*/
/*
void LinearAlgebra::multiplyMatrixA(vector<double>& x, vector<double>& b) {


//      [ H            -omega * I ] [X_real]  = [ b ]
//      [ omega * I     H         ] [X_imag]

//      double               mOmega;
//      vector<OuterProduct> mHessianMatrix;

if(just_once){
  cout << "Omega = " << mOmega << endl;
  cout << "x = [" << x[0];
  for(int cnt = 1; cnt < mRealDimension; ++cnt)
  cout << "; " << x[cnt];
  cout << "];" << endl;
}

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

        if(just_once){




          cout << "B = " << op.mScalar << " * ["
               << op.mDeltaPosition[0] * op.mDeltaPosition[0] << ", "
               << op.mDeltaPosition[0] * op.mDeltaPosition[1] << ", "
               << op.mDeltaPosition[0] * op.mDeltaPosition[2] << ";"
               << op.mDeltaPosition[1] * op.mDeltaPosition[0] << ", "
               << op.mDeltaPosition[1] * op.mDeltaPosition[1] << ", "
               << op.mDeltaPosition[1] * op.mDeltaPosition[2] << ";"
               << op.mDeltaPosition[2] * op.mDeltaPosition[0] << ", "
               << op.mDeltaPosition[2] * op.mDeltaPosition[1] << ", "
               << op.mDeltaPosition[2] * op.mDeltaPosition[2] << "];" << endl;

          cout << "A("     << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
               << ") = A(" << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
               << ") + B;" << endl;
          cout << "A("     << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
               << ") = A(" << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
               << ") - B;" << endl;
          cout << "A("     << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
               << ") = A(" << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex1 + 1 << ":" << op.mNodeIndex1 + 3
               << ") - B;" << endl;
          cout << "A("     << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
               << ") = A(" << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3 << "," << op.mNodeIndex2 + 1 << ":" << op.mNodeIndex2 + 3
               << ") + B;" << endl;

        }


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

    if(just_once){
      cout << "b = A * x" << endl;
    just_once = false;
    for(int cnt = 0; cnt < 2*mRealDimension; ++cnt)
    cout << "b[" << cnt << "] = " << b[cnt] << endl;
  }
}
*/
/*
void LinearAlgebra::removeProjIntoZeroEigen( vector<double>& u ) {
    double u_dot_xTran = 0.0;
    double u_dot_yTran = 0.0;
    double u_dot_zTran = 0.0;
    double u_dot_xyRot = 0.0;
    double u_dot_xzRot = 0.0;
    double u_dot_yzRot = 0.0;

    int index = 0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
        u_dot_xTran += u[index];
        u_dot_xyRot += u[index] * mRotationXY[index];
        u_dot_xzRot += u[index] * mRotationXZ[index];
        u_dot_yzRot += u[index] * mRotationYZ[index];

        ++index;
        u_dot_yTran += u[index];
        u_dot_xyRot += u[index] * mRotationXY[index];
        u_dot_xzRot += u[index] * mRotationXZ[index];
        u_dot_yzRot += u[index] * mRotationYZ[index];

        ++index;
        u_dot_zTran += u[index];
        u_dot_xzRot += u[index] * mRotationXZ[index];
        u_dot_yzRot += u[index] * mRotationYZ[index];

        ++index;
    }

    u_dot_xTran *= mNormalizedTranslation;
    u_dot_yTran *= mNormalizedTranslation;
    u_dot_zTran *= mNormalizedTranslation;

    index = 0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
        u[index] -= (u_dot_xTran + u_dot_xyRot * mRotationXY[index] + u_dot_xzRot * mRotationXZ[index] + u_dot_yzRot * mRotationYZ[index]);

        ++index;
        u[index] -= (u_dot_yTran + u_dot_xyRot * mRotationXY[index] + u_dot_xzRot * mRotationXZ[index] + u_dot_yzRot * mRotationYZ[index]);

        ++index;
        u[index] -= (u_dot_zTran                                    + u_dot_xzRot * mRotationXZ[index] + u_dot_yzRot * mRotationYZ[index]);

        ++index;
    }
}


void LinearAlgebra::CreateRotationalEigenVectors( vector<Node>& nodes ) {

    vector<double> xCoordinates;
    vector<double> yCoordinates;
    vector<double> zCoordinates;

    xCoordinates.reserve( mNumberNodes );
    yCoordinates.reserve( mNumberNodes );
    zCoordinates.reserve( mNumberNodes );

    double mean_x = 0.0;
    double mean_y = 0.0;
    double mean_z = 0.0;
    for(auto node : nodes) {
        xCoordinates.push_back( node.GetX() );
        yCoordinates.push_back( node.GetY() );
        zCoordinates.push_back( node.GetZ() );

        mean_x += node.GetX();
        mean_y += node.GetY();
        mean_z += node.GetZ();
    }

    mean_x *= mNormalizedTranslation;
    mean_y *= mNormalizedTranslation;
    mean_z *= mNormalizedTranslation;

    double xSq = 0.0;
    double ySq = 0.0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
      xCoordinates[cnt] -= mean_x;
      yCoordinates[cnt] -= mean_y;
      zCoordinates[cnt] -= mean_z;

      xSq += xCoordinates[cnt] * xCoordinates[cnt];
      ySq += yCoordinates[cnt] * yCoordinates[cnt];
    }

    double normalization_term = 1.0 / sqrt( xSq + ySq );
    double dot_xy_xz = 0.0;
    double dot_xy_yz = 0.0;
    int    index = 0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
        // Orthogonal to translation eigen vectors and normalized
        mRotationXY.push_back(  yCoordinates[cnt] * normalization_term );
        mRotationXY.push_back( -xCoordinates[cnt] * normalization_term );
        mRotationXY.push_back(  0.0 );

        // Following two eigen vectors are orthogonal to translation eigen vectors but not normalized
        mRotationXZ.push_back(  zCoordinates[cnt] );
        mRotationXZ.push_back(  0.0 );
        mRotationXZ.push_back( -xCoordinates[cnt] );
        index = 3*cnt;
        dot_xy_xz += mRotationXY[index] * mRotationXZ[index];

        mRotationYZ.push_back(  0.0 );
        mRotationYZ.push_back(  zCoordinates[cnt] );
        mRotationYZ.push_back( -yCoordinates[cnt] );
        ++index;
        dot_xy_yz += mRotationXY[index] * mRotationYZ[index];
    }

    normalization_term = xSq; // sym of squares for mRotationXZ z-component
    double dot_xz_yz   = 0.0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
        index = 3*cnt;
        mRotationXZ[index]  -=  dot_xy_xz * mRotationXY[index];
        mRotationYZ[index]   = -dot_xy_yz * mRotationXY[index];
        normalization_term  +=  mRotationXZ[index] * mRotationXZ[index];
        dot_xz_yz           +=  mRotationXZ[index] * mRotationYZ[index];

        ++index;
        mRotationXZ[index]  = -dot_xy_xz * mRotationXY[index];
        mRotationYZ[index] -=  dot_xy_yz * mRotationXY[index];

        normalization_term  +=  mRotationXZ[index] * mRotationXZ[index];
        dot_xz_yz           +=  mRotationXZ[index] * mRotationYZ[index];

        ++index;
        dot_xz_yz           +=  mRotationXZ[index] * mRotationYZ[index];
    }

    // normalize mRotationXZ and orthogonalize mRotationYZ
    normalization_term  = 1.0 / sqrt( normalization_term );
    dot_xz_yz          *= normalization_term;
    double sum_of_sq    = 0.0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
        index = 3*cnt;
        mRotationXZ[index] *= normalization_term;
        mRotationYZ[index] -= dot_xz_yz * mRotationXZ[index];
        sum_of_sq          += mRotationYZ[index] * mRotationYZ[index];

        ++index;
        mRotationXZ[index] *= normalization_term;
        mRotationYZ[index] -= dot_xz_yz * mRotationXZ[index];
        sum_of_sq          += mRotationYZ[index] * mRotationYZ[index];

        ++index;
        mRotationXZ[index] *= normalization_term;
        mRotationYZ[index] -= dot_xz_yz * mRotationXZ[index];
        sum_of_sq          += mRotationYZ[index] * mRotationYZ[index];
    }

    normalization_term  = 1.0 / sqrt( sum_of_sq );
    index = 0;
    for(int cnt = 0; cnt < mNumberNodes; ++cnt) {
        mRotationYZ[index++] *= normalization_term;
        mRotationYZ[index++] *= normalization_term;
        mRotationYZ[index++] *= normalization_term;
    }
}
*/
