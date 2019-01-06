#ifndef ____ENM__
#define ____ENM__

#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <vector>

#include "network_model/Connection.h"
#include "network_model/Node.h"
#include "network_model/HessianMatrix.h"
#include "network_model/ReferencePoint.h"
#include "util/CommonType.h"
#include "util/Common.h"

using namespace std;

class ElasticNetworkModel{
public:
    ElasticNetworkModel();
    ~ElasticNetworkModel();

    void SetMaxSpringLength( double m ) { mMaxSpringLength = m; mMaxRangeSq = m*m; }
    void SetNodeType(string);

    void AllocateNodes(int s) { mNodes.reserve(s); }

    int AddNode(
      int amino_id, string amino_name,
      vector<int>, vector<string>, vector<double>, vector<double>, vector<double>, vector<double>);

    void CreateBackboneConnection(int, int);
    void IdentifyContacts();

    void ConstructLinearResponse();
    void SetConstantVector(vector<double>& b);

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

    // Used internally as part of node construction
    double mWeight;
    double mPosX;
    double mPosY;
    double mPosZ;
    double mEstQ;

    vector<ReferencePoint> mReferencePoints;

    HessianMatrix mHessianMatrix;

    Common mParameter;
};

#endif
