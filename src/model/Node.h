#ifndef ____Node__
#define ____Node__

#include <iostream>

using namespace std;

class Node{
public:
    Node(int id, double x, double y, double z);
    ~Node();
    
    void SetPosition(double x, double y, double z);
    void SetIndex   (int    x, int    y, int    z);
    
    double GetX() { return mPos[0]; }
    double GetY() { return mPos[1]; }
    double GetZ() { return mPos[2]; }
    int    GetID() { return mNodeID; }
    
    
private:
    int    mNodeID;
    int    mVoxelIndex[3];
    double mPos[3];
    
};

#endif
