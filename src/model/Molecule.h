#ifndef ____Molecule__
#define ____Molecule__

#include <iostream>
#include <vector>

#include "../../src/model/Node.h"

using namespace std;

class Molecule{
public:
    Molecule();
    ~Molecule();
    
    void AllocateNodes(int s) { mNodes.reserve(s); }
    void AddNode(double x, double y, double z);
    void IdentifyContacts( double cutoff );
    
    
    void PrintLimits();
    
private:
    vector<Node> mNodes;
    
    double mMinX;
    double mMaxX;
    double mMinY;
    double mMaxY;
    double mMinZ;
    double mMaxZ;
    /*
     double mMinX =  numeric_limits<double>::max();
     double mMaxX = -numeric_limits<double>::max();//numeric_limits<double>::lowest();
     double mMinY =  numeric_limits<double>::max();
     double mMaxY = -numeric_limits<double>::max();//numeric_limits<double>::lowest();
     double mMinZ =  numeric_limits<double>::max();
     double mMaxZ = -numeric_limits<double>::max();//numeric_limits<double>::lowest();
     */
};

#endif
