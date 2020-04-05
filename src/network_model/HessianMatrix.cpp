#include <fstream>

#include "network_model/HessianMatrix.h"

HessianMatrix::HessianMatrix(){
  mNullSpacePtr = nullptr;
}

HessianMatrix::~HessianMatrix(){
}

void HessianMatrix::SetConnection(int node_1, int node_2, vector<double> dp, double k) {
  mHessianMatrix.push_back( {3*node_1, 3*node_2, {dp[0], dp[1], dp[2]}, k} );
}

void HessianMatrix::MultiplyMatrix(vector<double>& x, vector<double>& b) {

    for(int cnt = 0; cnt < static_cast<int>(x.size()); ++cnt) {
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

  //  if(mNullSpacePtr != nullptr) {
  //    mNullSpacePtr->RemoveProjection( b );
  //  }
}


void HessianMatrix::PrintHessian(){

  int max_index = 0;
  for(OuterProduct op : mHessianMatrix) {
    if( max_index <  op.mNodeIndex1 )
      max_index =  op.mNodeIndex1;

    if( max_index <  op.mNodeIndex2 )
      max_index =  op.mNodeIndex2;
  }
  max_index += 3;
  vector<double> hessian_matrix(max_index*max_index, 0.0);

  for(OuterProduct op : mHessianMatrix) {
    double xx = op.mScalar * op.mDeltaPosition[0] * op.mDeltaPosition[0];
    double xy = op.mScalar * op.mDeltaPosition[0] * op.mDeltaPosition[1];
    double xz = op.mScalar * op.mDeltaPosition[0] * op.mDeltaPosition[2];
    double yy = op.mScalar * op.mDeltaPosition[1] * op.mDeltaPosition[1];
    double yz = op.mScalar * op.mDeltaPosition[1] * op.mDeltaPosition[2];
    double zz = op.mScalar * op.mDeltaPosition[2] * op.mDeltaPosition[2];

    hessian_matrix[(op.mNodeIndex1 + 0) * max_index + op.mNodeIndex1 + 0] += xx;
    hessian_matrix[(op.mNodeIndex1 + 0) * max_index + op.mNodeIndex1 + 1] += xy;
    hessian_matrix[(op.mNodeIndex1 + 0) * max_index + op.mNodeIndex1 + 2] += xz;
    hessian_matrix[(op.mNodeIndex1 + 1) * max_index + op.mNodeIndex1 + 0] += xy;
    hessian_matrix[(op.mNodeIndex1 + 1) * max_index + op.mNodeIndex1 + 1] += yy;
    hessian_matrix[(op.mNodeIndex1 + 1) * max_index + op.mNodeIndex1 + 2] += yz;
    hessian_matrix[(op.mNodeIndex1 + 2) * max_index + op.mNodeIndex1 + 0] += xz;
    hessian_matrix[(op.mNodeIndex1 + 2) * max_index + op.mNodeIndex1 + 1] += yz;
    hessian_matrix[(op.mNodeIndex1 + 2) * max_index + op.mNodeIndex1 + 2] += zz;

    hessian_matrix[(op.mNodeIndex2 + 0) * max_index + op.mNodeIndex1 + 0] -= xx;
    hessian_matrix[(op.mNodeIndex2 + 0) * max_index + op.mNodeIndex1 + 1] -= xy;
    hessian_matrix[(op.mNodeIndex2 + 0) * max_index + op.mNodeIndex1 + 2] -= xz;
    hessian_matrix[(op.mNodeIndex2 + 1) * max_index + op.mNodeIndex1 + 0] -= xy;
    hessian_matrix[(op.mNodeIndex2 + 1) * max_index + op.mNodeIndex1 + 1] -= yy;
    hessian_matrix[(op.mNodeIndex2 + 1) * max_index + op.mNodeIndex1 + 2] -= yz;
    hessian_matrix[(op.mNodeIndex2 + 2) * max_index + op.mNodeIndex1 + 0] -= xz;
    hessian_matrix[(op.mNodeIndex2 + 2) * max_index + op.mNodeIndex1 + 1] -= yz;
    hessian_matrix[(op.mNodeIndex2 + 2) * max_index + op.mNodeIndex1 + 2] -= zz;

    hessian_matrix[(op.mNodeIndex1 + 0) * max_index + op.mNodeIndex2 + 0] -= xx;
    hessian_matrix[(op.mNodeIndex1 + 0) * max_index + op.mNodeIndex2 + 1] -= xy;
    hessian_matrix[(op.mNodeIndex1 + 0) * max_index + op.mNodeIndex2 + 2] -= xz;
    hessian_matrix[(op.mNodeIndex1 + 1) * max_index + op.mNodeIndex2 + 0] -= xy;
    hessian_matrix[(op.mNodeIndex1 + 1) * max_index + op.mNodeIndex2 + 1] -= yy;
    hessian_matrix[(op.mNodeIndex1 + 1) * max_index + op.mNodeIndex2 + 2] -= yz;
    hessian_matrix[(op.mNodeIndex1 + 2) * max_index + op.mNodeIndex2 + 0] -= xz;
    hessian_matrix[(op.mNodeIndex1 + 2) * max_index + op.mNodeIndex2 + 1] -= yz;
    hessian_matrix[(op.mNodeIndex1 + 2) * max_index + op.mNodeIndex2 + 2] -= zz;

    hessian_matrix[(op.mNodeIndex2 + 0) * max_index + op.mNodeIndex2 + 0] += xx;
    hessian_matrix[(op.mNodeIndex2 + 0) * max_index + op.mNodeIndex2 + 1] += xy;
    hessian_matrix[(op.mNodeIndex2 + 0) * max_index + op.mNodeIndex2 + 2] += xz;
    hessian_matrix[(op.mNodeIndex2 + 1) * max_index + op.mNodeIndex2 + 0] += xy;
    hessian_matrix[(op.mNodeIndex2 + 1) * max_index + op.mNodeIndex2 + 1] += yy;
    hessian_matrix[(op.mNodeIndex2 + 1) * max_index + op.mNodeIndex2 + 2] += yz;
    hessian_matrix[(op.mNodeIndex2 + 2) * max_index + op.mNodeIndex2 + 0] += xz;
    hessian_matrix[(op.mNodeIndex2 + 2) * max_index + op.mNodeIndex2 + 1] += yz;
    hessian_matrix[(op.mNodeIndex2 + 2) * max_index + op.mNodeIndex2 + 2] += zz;
  }

  fstream fs;
  fs.open ("Hessian.txt", std::fstream::out);

  for(int r = 0; r < max_index; ++r){
    for(int c = 0; c < max_index; ++c){
      fs << r << " " << c << " " << hessian_matrix[r * max_index + c] << endl;
    }
  }

  fs.close();
}



void HessianMatrix::PrintMatrixA(){

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
