#ifndef ____Node__
#define ____Node__

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include "util/CommonType.h"
#include "util/Common.h"

using namespace std;

class Node{
public:
    Node(Common&,
         int    id,
         int    amino_id,
         string amino_name,
         double x,
         double y,
         double z,
         double q,
         vector<int>,
         vector<string>,
         vector<double>,
         vector<double>,
         vector<double>,
         vector<double>);
    ~Node();

    void SetPosition(double x, double y, double z);

    double GetX() { return mPos[0]; }
    double GetY() { return mPos[1]; }
    double GetZ() { return mPos[2]; }
    double GetQ() { return mCharge; }
    int    GetID() { return mNodeID; }

    vector<double> GetSeparationVector(Node&);

    void   AddConnection(CommonType::ConnectionType, int);
    bool   IsConnected(int, CommonType::ConnectionType);

    vector<double> GetElectricPotential( vector<double>& probe_loc );


    void Print();

private:
    int                  mNodeID;
    double               mPos[3];
    double               mCharge;

    int                  mResidualId;
    string               mResidualName;

    vector<int>          mAtomId;
    vector<string>       mAtomName;
    vector<double>       mAtomPosX;
    vector<double>       mAtomPosY;
    vector<double>       mAtomPosZ;
    vector<double>       mAtomCharge;

    // Maps a connection type to a set of node IDs this node is connected to
    map< CommonType::ConnectionType, set<int> > mConnections;
    // set<int> nodeIds = mConnections[CommonEnum::ConnectionType::SPRING_LEVEL_1];

    Common& mParameter;
};

#endif
