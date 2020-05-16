#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>
#include "Node.h"

using namespace std;

class Kernel {

public:

  Kernel( vector<Node>& node_position ){
    mSize = node_position.size();

    mTranslation = 1. / static_cast<double>(mSize);
    mRotationXY  = new double[ 2 * node_position.size() ]();
    mRotationXZ  = new double[ 3 * node_position.size() ]();
    mRotationYZ  = new double[ 3 * node_position.size() ]();

    double* node_pos;
    double  dot_product[6] = {0. , 0., 0., 0. , 0., 0.};
    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      node_pos = node_position[cnt].GetPosition();

      dot_product[0] += node_pos[0];
      dot_product[1] += node_pos[1];
      dot_product[2] += node_pos[2];

      mRotationXY[2*cnt  ] =  node_pos[1];
      mRotationXY[2*cnt+1] = -node_pos[0];

      mRotationXZ[3*cnt  ] =  node_pos[2];
      mRotationXZ[3*cnt+2] = -node_pos[0];

      mRotationYZ[3*cnt+1] =  node_pos[2];
      mRotationYZ[3*cnt+2] = -node_pos[1];
    }

    dot_product[0] *= mTranslation;
    dot_product[1] *= mTranslation;
    dot_product[2] *= mTranslation;

    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      mRotationXY[2*cnt  ] -= dot_product[1];
      dot_product[3] += mRotationXY[2*cnt  ] * mRotationXY[2*cnt  ];
      mRotationXY[2*cnt+1] += dot_product[0];
      dot_product[3] += mRotationXY[2*cnt+1] * mRotationXY[2*cnt+1];

      mRotationXZ[3*cnt  ] -= dot_product[2];
      dot_product[4] += mRotationXY[2*cnt  ] * mRotationXZ[3*cnt  ];
      mRotationXZ[3*cnt+2] += dot_product[0];

      mRotationYZ[3*cnt+1] -= dot_product[2];
      dot_product[5] += mRotationXY[2*cnt+1] * mRotationYZ[3*cnt+1];
      mRotationYZ[3*cnt+2] += dot_product[1];
    }

    dot_product[3] = sqrt( dot_product[3] );
    dot_product[4] /= dot_product[3];
    dot_product[5] /= dot_product[3];

    dot_product[0] = 0.;
    dot_product[1] = 0.;
    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      mRotationXY[2*cnt  ] /= dot_product[3];
      mRotationXY[2*cnt+1] /= dot_product[3];

      mRotationXZ[3*cnt  ] -= mRotationXY[2*cnt  ] * dot_product[4];
      dot_product[0] += mRotationXZ[3*cnt  ] * mRotationXZ[3*cnt  ];
      mRotationXZ[3*cnt+1] -= mRotationXY[2*cnt+1] * dot_product[4];
      dot_product[0] += mRotationXZ[3*cnt+1] * mRotationXZ[3*cnt+1] + mRotationXZ[3*cnt+2] * mRotationXZ[3*cnt+2];

      mRotationYZ[3*cnt  ] -= mRotationXY[2*cnt  ] * dot_product[5];
      dot_product[1] += mRotationXZ[3*cnt  ] * mRotationYZ[3*cnt  ];
      mRotationYZ[3*cnt+1] -= mRotationXY[2*cnt+1] * dot_product[5];
      dot_product[1] += mRotationXZ[3*cnt+1] * mRotationYZ[3*cnt+1] + mRotationXZ[3*cnt+2] * mRotationYZ[3*cnt+2];
    }

    dot_product[0] = sqrt( dot_product[0] );
    dot_product[1] /= dot_product[0];

    dot_product[2] = 0.;
    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      mRotationXZ[3*cnt  ] /= dot_product[0];
      mRotationXZ[3*cnt+1] /= dot_product[0];
      mRotationXZ[3*cnt+2] /= dot_product[0];

      mRotationYZ[3*cnt  ] -= mRotationXZ[3*cnt  ] * dot_product[1];
      dot_product[2] += mRotationYZ[3*cnt  ] * mRotationYZ[3*cnt  ];
      mRotationYZ[3*cnt+1] -= mRotationXZ[3*cnt+1] * dot_product[1];
      dot_product[2] += mRotationYZ[3*cnt+1] * mRotationYZ[3*cnt+1];
      mRotationYZ[3*cnt+2] -= mRotationXZ[3*cnt+2] * dot_product[1];
      dot_product[2] += mRotationYZ[3*cnt+2] * mRotationYZ[3*cnt+2];
    }

    dot_product[2] = sqrt( dot_product[2] );

    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      mRotationYZ[3*cnt  ] /= dot_product[2];
      mRotationYZ[3*cnt+1] /= dot_product[2];
      mRotationYZ[3*cnt+2] /= dot_product[2];
    }
  }
  ~Kernel(){
    delete[] mRotationXY;
    delete[] mRotationXZ;
    delete[] mRotationYZ;
  }

  void RemoveProjection( double* v ){
    double  dot_product[6] = {0. , 0., 0., 0. , 0., 0.};
    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      dot_product[0] += v[3*cnt  ];
      dot_product[1] += v[3*cnt+1];
      dot_product[2] += v[3*cnt+2];

      dot_product[3] += v[3*cnt] * mRotationXY[2*cnt] + v[3*cnt+1] * mRotationXY[2*cnt+1];
      dot_product[4] += v[3*cnt] * mRotationXZ[3*cnt] + v[3*cnt+1] * mRotationXZ[3*cnt+1] + v[3*cnt+2] * mRotationXZ[3*cnt+2];
      dot_product[5] += v[3*cnt] * mRotationYZ[3*cnt] + v[3*cnt+1] * mRotationYZ[3*cnt+1] + v[3*cnt+2] * mRotationYZ[3*cnt+2];
    }

    dot_product[0] *= mTranslation;
    dot_product[1] *= mTranslation;
    dot_product[2] *= mTranslation;

    for( size_t cnt = 0; cnt < mSize; ++cnt ){
      v[3*cnt  ] -= (dot_product[0] + dot_product[3]*mRotationXY[2*cnt  ] + dot_product[4]*mRotationXZ[3*cnt  ]  + dot_product[5]*mRotationYZ[3*cnt  ] );
      v[3*cnt+1] -= (dot_product[1] + dot_product[3]*mRotationXY[2*cnt+1] + dot_product[4]*mRotationXZ[3*cnt+1]  + dot_product[5]*mRotationYZ[3*cnt+1] );
      v[3*cnt+2] -= (dot_product[2]                                       + dot_product[4]*mRotationXZ[3*cnt+2]  + dot_product[5]*mRotationYZ[3*cnt+2] );
    }
  }

private:

  size_t  mSize = 0;

  double  mTranslation = 1.0;
  double* mRotationXY;
  double* mRotationXZ;
  double* mRotationYZ;

};

#endif
