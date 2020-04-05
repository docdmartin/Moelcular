#ifndef ____HessianMatrix__
#define ____HessianMatrix__

#include <iostream>
#include <vector>

#include "network_model/NullSpace.h"

using namespace std;

class HessianMatrix{
public:
    HessianMatrix();
    ~HessianMatrix();

    void AddNullSpace( NullSpace* np ) { mNullSpacePtr = np; }
    void RemoveProjection( vector<double>& a ) { if( mNullSpacePtr != nullptr ){ mNullSpacePtr->RemoveProjection( a ); } }

    void SetConnection(int, int, vector<double>, double);
    void MultiplyMatrix(vector<double>&, vector<double>&);

    void PrintMatrixA();
    void PrintHessian();

private:
    struct OuterProduct {
        int    mNodeIndex1;
        int    mNodeIndex2;

        /* position[mNodeIndex1] - position[mNodeIndex2] */
        double mDeltaPosition[3];

        /* spring constant for this specific node pair */
        double mScalar;
    };

    vector<OuterProduct> mHessianMatrix;

    NullSpace* mNullSpacePtr;
};

#endif
