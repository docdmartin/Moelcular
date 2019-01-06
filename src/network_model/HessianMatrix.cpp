#include "network_model/HessianMatrix.h"

HessianMatrix::HessianMatrix(){
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
