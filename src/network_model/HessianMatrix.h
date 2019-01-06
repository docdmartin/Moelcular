#ifndef ____HessianMatrix__
#define ____HessianMatrix__

#include <iostream>
#include <vector>

using namespace std;

class HessianMatrix{
public:
    HessianMatrix();
    ~HessianMatrix();

    void SetConnection(int, int, vector<double>, double);
    void MultiplyMatrix(vector<double>&, vector<double>&);

    void PrintMatrixA();

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

};

#endif
