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

    void Print();

private:
    Node& mNode1;
    Node& mNode2;

    CommonType::ConnectionType mConnectionType;
    vector<double> mSeparation;
    double         mSquaredDistance;

    Common& mParameter;
};

#endif
