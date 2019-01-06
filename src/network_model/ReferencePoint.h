#ifndef ____REFPNT__
#define ____REFPNT__

#include <iostream>
#include <vector>

#include "network_model/Node.h"
#include "network_model/HessianMatrix.h"

using namespace std;

class ReferencePoint{
public:
    ReferencePoint(double x, double y, double z);
    ~ReferencePoint();

    void BuildPotential(vector<Node> &, HessianMatrix &);

    pair<double, double> SingleModeFrequencyResponse(double omega);
    pair<double, double> DualModeFrequencyResponse(double omega1, double omega2, double a);

    void Print();

private:
    double innerProduct( vector<double>&, vector<double>& );
    void   grahamSchmit( vector<double>&, vector<double>&, double );
    double rescale     ( vector<double>& );
    void   partialSolver(vector<double>&, vector<double>&, double*, double*);

    vector<double> mPosition;
    double         mBTb;
    double         mSigma;
    vector<double> mLambda;
};

#endif
