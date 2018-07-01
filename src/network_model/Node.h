#ifndef ____Node__
#define ____Node__

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include "util/Common.h"

using namespace std;

class Node{
public:
    Node(int id, double x, double y, double z);
    ~Node();

    void SetPosition(double x, double y, double z);

    double GetX() { return mPos[0]; }
    double GetY() { return mPos[1]; }
    double GetZ() { return mPos[2]; }
    int    GetID() { return mNodeID; }

    vector<double> GetSeparationVector(Node&);

    void   AddConnection(CommonEnum::ConnectionType, int);
    bool   IsConnected(int, CommonEnum::ConnectionType);


    void Print();

private:
    int    mNodeID;
    CommonEnum::NodeType mNodeType = CommonEnum::NodeType::C_ALPHA;
    double mPos[3];

    map< CommonEnum::ConnectionType, set<int> > mConnections;
};

#endif
