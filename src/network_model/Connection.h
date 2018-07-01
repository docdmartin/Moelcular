#ifndef ____Connection__
#define ____Connection__

#include <iostream>
#include <vector>

#include "network_model/Connection.h"
#include "network_model/Node.h"
#include "util/Common.h"

using namespace std;

class Connection{
public:
    Connection(Node&, Node&, CommonEnum::ConnectionType);
    ~Connection();

    void Print();

private:
    Node& mNode1;
    Node& mNode2;

    CommonEnum::ConnectionType mConnectionType;
    vector<double> mSeparation;
    double         mSquaredDistance;
};

#endif
