#ifndef ____Connection__
#define ____Connection__

#include <iostream>
#include <vector>

#include "network_model/Node.h"

using namespace std;

class Connection{
public:
    Connection(Common&, Node&, Node&, CommonType::ConnectionType);
    ~Connection();

    int GetNodeId1() { return mNode1.GetID(); }
    int GetNodeId2() { return mNode2.GetID(); }

    vector<double> GetSeparation() { return mSeparation; }

    double GetSpringConstant() { return mSpringConstant; }

    void Print();

private:
    Node& mNode1;
    Node& mNode2;

    CommonType::ConnectionType mConnectionType;
    vector<double> mSeparation;
    double         mSquaredDistance;
    double         mSpringConstant;

    Common& mParameter;
};

#endif
