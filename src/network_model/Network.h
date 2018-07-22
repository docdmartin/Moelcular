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

    int AddNode(
      int amino_id, string amino_name,
      vector<int>, vector<string>, vector<double>, vector<double>, vector<double>, vector<double>);

    void CreateBackboneConnection(int, int);
    void IdentifyContacts( double cutoff );


    void Print();

private:
    void addConnection(CommonEnum::ConnectionType, Node&, Node&);
    bool calculateNodeProperties(vector<string>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

    CommonEnum::NodeType                                mEstimateNodeType;
    vector<Node>                                        mNodes;
    map<CommonEnum::ConnectionType, vector<Connection>> mConnections;

    // Used internally as part of node construction
    double mWeight;
    double mPosX;
    double mPosY;
    double mPosZ;
    double mEstQ;
};

#endif
