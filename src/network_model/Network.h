#ifndef ____Network__
#define ____Network__

#include <iostream>
#include <vector>

#include "network_model/Node.h"
#include "network_model/Connection.h"
#include "util/Common.h"

using namespace std;

class Network{
public:
    Network();
    ~Network();

    void AllocateNodes(int s) { mNodes.reserve(s); }

    void AddNode(double x, double y, double z);
    void AddConnection(CommonEnum::ConnectionType, Node&, Node&);

    void IdentifyContacts( double cutoff );


    void Print();

private:
    vector<Node> mNodes;
    map<CommonEnum::ConnectionType, vector<Connection>> mConnections;
};

#endif
