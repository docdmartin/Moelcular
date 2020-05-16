#ifndef MATRIX_3X3_H
#define MATRIX_3X3_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "Common.h"

using namespace std;

class Matrix3x3 {

public:

  Matrix3x3( size_t i1, size_t i2, double* p1, double* p2, double k )
  {
    mIndex[0] = 3*i1;
    mIndex[1] = 3*i2;

    double x = p1[0] - p2[0];
    double y = p1[1] - p2[1];
    double z = p1[2] - p2[2];
    double mag = std::sqrt( x*x + y*y + z*z );
    mRange = mag;
    if( mag == 0.0 )
      throw string("Neighboring points have zero magnitude");
    mUnitVector[0] = x / mag;
    mUnitVector[1] = y / mag;
    mUnitVector[2] = z / mag;
    mDotVector[0] = k*mUnitVector[0];
    mDotVector[1] = k*mUnitVector[1];
    mDotVector[2] = k*mUnitVector[2];
  }
  ~Matrix3x3(){}

  void Multiply( double* b, double* x )
  {
    double dot_prod =
      mDotVector[0] * ( x[mIndex[1]  ] - x[mIndex[0]  ] ) +
      mDotVector[1] * ( x[mIndex[1]+1] - x[mIndex[0]+1] ) +
      mDotVector[2] * ( x[mIndex[1]+2] - x[mIndex[0]+2] );

    double dx = dot_prod * mUnitVector[0];
    double dy = dot_prod * mUnitVector[1];
    double dz = dot_prod * mUnitVector[2];

    b[mIndex[0]  ] -= dx;
    b[mIndex[0]+1] -= dy;
    b[mIndex[0]+2] -= dz;

    b[mIndex[1]  ] += dx;
    b[mIndex[1]+1] += dy;
    b[mIndex[1]+2] += dz;
  }

  void Print(){
    if( mIndex[0] != 0 )
      return;

    cout << "Index: " << mIndex[0] << ", " << mIndex[1] << ", range: " << mRange << ", upper: "
    << mUnitVector[0] * mDotVector[0] << ", " << mUnitVector[0] * mDotVector[1] << ", " << mUnitVector[0] * mDotVector[2] << ", "
    << mUnitVector[1] * mDotVector[1] << ", " << mUnitVector[1] * mDotVector[2] << ", "
    << mUnitVector[2] * mDotVector[2] << endl;
  }

private:

  double      mDotVector[3];
  double      mUnitVector[3];
  size_t      mIndex[2];

  double      mRange;

};

#endif
