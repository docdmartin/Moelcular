#ifndef ____NULLSPACE__
#define ____NULLSPACE__

#include <iostream>
#include <vector>

#include "network_model/Node.h"

using namespace std;

class NullSpace{
public:
    NullSpace(){}
    ~NullSpace(){}

    void CreateEigenVectors(vector<Node>& nodes);
    void RemoveProjection( vector<double>& );

private:
    vector<double> mNullEigenvectors[6];
    int mDimension;
};

#endif
