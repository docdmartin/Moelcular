#ifndef ____LinearAlgebra__
#define ____LinearAlgebra__

#include <iostream>
#include <vector>

#include "network_model/Node.h"

using namespace std;

class LinearAlgebra{
public:
    LinearAlgebra();
    ~LinearAlgebra();

    void SetDimension(int num_of_nodes);

    void PreProcess();
    void PerformCG();
    void GetPotential(double* real_x, double* imag_x);

    void SetOmega(double d) { mOmega = d; }
    void SetConnection(int, int, vector<double>, double);
    void SetConstantVector(vector<double>&);
//    void CreateRotationalEigenVectors( vector<Node>& );

private:
    struct OuterProduct {
        int    mNodeIndex1;
        int    mNodeIndex2;

        /* position[mNodeIndex1] - position[mNodeIndex2] */
        double mDeltaPosition[3];

        /* spring constant for this specific node pair */
        double mScalar;
    };

    bool just_once = true;

    void multiplyLinearSystem(vector<double>&, vector<double>&);
    void multiplyMatrixA(vector<double>&, vector<double>&);
//    void removeProjIntoZeroEigen( vector<double>& );

    void printMatrixA();
    void printMatrixB();


    int    mNumberNodes;
    int    mRealDimension;
    size_t mSize;

    double mToleranceSq;
    double mSolverConvergeThreshold;

/*
    [ H            -omega * I ] [X_real]  = [ b ]
    [ omega * I     H         ] [X_imag]    [ 0 ]
*/
    double               mOmega;
    vector<OuterProduct> mHessianMatrix;
    vector<double      > mSolution;
    vector<double      > mConstantTerm;

    // All used in CG Solver
    vector<double> B;
    vector<double> r;
    vector<double> z;
    vector<double> p;
    vector<double> w;


    // translation eigenvector
    double mNormalizedTranslation;
    vector<double> mRotationXY;
    vector<double> mRotationXZ;
    vector<double> mRotationYZ;
};

#endif
