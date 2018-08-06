#ifndef ____Network__
#define ____Network__

#define X_INDEX   0
#define Y_INDEX   1
#define Z_INDEX   2
#define MIN_INDEX 0
#define MAX_INDEX 1

#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <vector>

#include "network_model/Connection.h"

using namespace std;

class Network{
public:
    Network(Common&);
    ~Network();

    void SetMaxSpringLength( double m ) { mMaxSpringLength = m; mMaxRangeSq = m*m; }

    void AllocateNodes(int s) { mNodes.reserve(s); }

    int AddNode(
      int amino_id, string amino_name,
      vector<int>, vector<string>, vector<double>, vector<double>, vector<double>, vector<double>);

    void CreateBackboneConnection(int, int);
    void IdentifyContacts();

    void Print();

private:
    void addConnection(CommonType::ConnectionType, Node&, Node&);
    void testConnection(int, int );
    bool calculateNodeProperties(vector<string>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);
    void voxelizeNode();

    double                                              mMaxSpringLength;
    double                                              mMaxRangeSq;
    vector<Node>                                        mNodes;
    map<CommonType::ConnectionType, vector<Connection>> mConnections;
    map<string, CommonType::ElementType>                mKnownElements;

    map<int, map<int, map<int, vector<int> > > > mNodeVoxels;

    double mAABB[3][2] = {{numeric_limits<double>::max(), numeric_limits<double>::lowest()},
                          {numeric_limits<double>::max(), numeric_limits<double>::lowest()},
                          {numeric_limits<double>::max(), numeric_limits<double>::lowest()}};

    // Used internally as part of node construction
    double mWeight;
    double mPosX;
    double mPosY;
    double mPosZ;
    double mEstQ;

    Common& mParameter;
};

#endif
