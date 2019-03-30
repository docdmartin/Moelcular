#include <cmath>

#include "network_model/NullSpace.h"
#include "network_model/Node.h"

using namespace std;



void NullSpace::CreateEigenVectors(vector<Node>& nodes) {
  mDimension = 3 * static_cast<int>( nodes.size() );
  for(int cnt = 0; cnt < 6; ++cnt) {
    mNullEigenvectors[cnt].resize( mDimension );
  }

  double v = 1.0 / sqrt( static_cast<double>( nodes.size() ) );

  double proj_3_on_0 = 0.0;
  double proj_3_on_1 = 0.0;

  double proj_4_on_0 = 0.0;
  double proj_4_on_2 = 0.0;

  double proj_5_on_1 = 0.0;
  double proj_5_on_2 = 0.0;

  double tmp;

  int n = -1;
  for( Node& node : nodes ) {

    tmp = v * node.GetX();
    proj_3_on_1 -= tmp;
    proj_4_on_2 -= tmp;

    tmp = v * node.GetY();
    proj_3_on_0 += tmp;
    proj_5_on_2 -= tmp;

    tmp = v * node.GetZ();
    proj_4_on_0 += tmp;
    proj_5_on_1 += tmp;

    ++n;
    mNullEigenvectors[0][n] = v;
    mNullEigenvectors[1][n] = 0.0;
    mNullEigenvectors[2][n] = 0.0;
    mNullEigenvectors[3][n] = node.GetY();
    mNullEigenvectors[4][n] = node.GetZ();
    mNullEigenvectors[5][n] = 0.0;

    ++n;
    mNullEigenvectors[0][n] = 0.0;
    mNullEigenvectors[1][n] = v;
    mNullEigenvectors[2][n] = 0.0;
    mNullEigenvectors[3][n] = -node.GetX();
    mNullEigenvectors[4][n] = 0.0;
    mNullEigenvectors[5][n] = node.GetZ();

    ++n;
    mNullEigenvectors[0][n] = 0.0;
    mNullEigenvectors[1][n] = 0.0;
    mNullEigenvectors[2][n] = v;
    mNullEigenvectors[3][n] = 0.0;
    mNullEigenvectors[4][n] = -node.GetX();
    mNullEigenvectors[5][n] = -node.GetY();
  }

  for(n = 0; n < mDimension; ++n) {
    switch( n % 3 )
    {
      case 0:
        mNullEigenvectors[3][n] -= proj_3_on_0 * mNullEigenvectors[0][n];
        mNullEigenvectors[4][n] -= proj_4_on_0 * mNullEigenvectors[0][n];
        break;
      case 1:
        mNullEigenvectors[3][n] -= proj_3_on_1 * mNullEigenvectors[1][n];
        mNullEigenvectors[5][n] -= proj_5_on_1 * mNullEigenvectors[1][n];
        break;
      case 2:
        mNullEigenvectors[4][n] -= proj_4_on_2 * mNullEigenvectors[2][n];
        mNullEigenvectors[5][n] -= proj_5_on_2 * mNullEigenvectors[2][n];
        break;
    }
  }

  // orthonormal 3:
  double proj_3_on_3 = 0.0;
  double proj_4_on_3 = 0.0;
  double proj_5_on_3 = 0.0;
  for(n = 0; n < mDimension; ++n) {
    proj_3_on_3 += mNullEigenvectors[3][n] * mNullEigenvectors[3][n];
    proj_4_on_3 += mNullEigenvectors[4][n] * mNullEigenvectors[3][n];
    proj_5_on_3 += mNullEigenvectors[5][n] * mNullEigenvectors[3][n];
  }

  proj_3_on_3 = 1.0 / sqrt( proj_3_on_3 );
  proj_4_on_3 *= proj_3_on_3;
  proj_5_on_3 *= proj_3_on_3;
  for(n = 0; n < mDimension; ++n) {
    mNullEigenvectors[3][n] *= proj_3_on_3;
    mNullEigenvectors[4][n] -= proj_4_on_3 * mNullEigenvectors[3][n];
    mNullEigenvectors[5][n] -= proj_5_on_3 * mNullEigenvectors[3][n];
  }

  // orthonormal 4:
  double proj_4_on_4 = 0.0;
  double proj_5_on_4 = 0.0;
  for(n = 0; n < mDimension; ++n) {
    proj_4_on_4 += mNullEigenvectors[4][n] * mNullEigenvectors[4][n];
    proj_5_on_4 += mNullEigenvectors[5][n] * mNullEigenvectors[4][n];
  }

  proj_4_on_4 = 1.0 / sqrt( proj_4_on_4 );
  proj_5_on_4 *= proj_4_on_4;
  for(n = 0; n < mDimension; ++n) {
    mNullEigenvectors[4][n] *= proj_4_on_4;
    mNullEigenvectors[5][n] -= proj_5_on_4 * mNullEigenvectors[4][n];
  }

  // normalize 5:
  double proj_5_on_5 = 0.0;
  for(n = 0; n < mDimension; ++n) {
    proj_5_on_5 += mNullEigenvectors[5][n] * mNullEigenvectors[5][n];
  }

  proj_5_on_5 = 1.0 / sqrt( proj_5_on_5 );
  for(n = 0; n < mDimension; ++n) {
    mNullEigenvectors[5][n] *= proj_5_on_5;
  }
}

void NullSpace::RemoveProjection( vector<double>& a ) {
  double proj_a_on_0 = 0.0;
  double proj_a_on_1 = 0.0;
  double proj_a_on_2 = 0.0;
  double proj_a_on_3 = 0.0;
  double proj_a_on_4 = 0.0;
  double proj_a_on_5 = 0.0;

  for(int n = 0; n < mDimension; ++n) {
    proj_a_on_0 += a[n] * mNullEigenvectors[0][n];
    proj_a_on_1 += a[n] * mNullEigenvectors[1][n];
    proj_a_on_2 += a[n] * mNullEigenvectors[2][n];
    proj_a_on_3 += a[n] * mNullEigenvectors[3][n];
    proj_a_on_4 += a[n] * mNullEigenvectors[4][n];
    proj_a_on_5 += a[n] * mNullEigenvectors[5][n];
  }

  for(int n = 0; n < mDimension; ++n) {
    a[n] -= proj_a_on_0 * mNullEigenvectors[0][n];
    a[n] -= proj_a_on_1 * mNullEigenvectors[1][n];
    a[n] -= proj_a_on_2 * mNullEigenvectors[2][n];
    a[n] -= proj_a_on_3 * mNullEigenvectors[3][n];
    a[n] -= proj_a_on_4 * mNullEigenvectors[4][n];
    a[n] -= proj_a_on_5 * mNullEigenvectors[5][n];
  }
}
