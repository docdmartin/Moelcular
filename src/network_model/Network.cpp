#include <map>
#include <limits>
#include <iostream>
#include <cmath>

#include "network_model/Network.h"

using namespace std;

struct coord {
    int ix, iy, iz;

    bool operator==(const coord &o) const {
        return (ix == o.ix && iy == o.iy && iz == o.iz);
    }

    bool operator<(const coord &o) const {
        return (ix < o.ix || (ix == o.ix && iy < o.iy) || (ix == o.ix && iy == o.iy && iz < o.iz));
    }
};



Network::Network() {
  mEstimateNodeType = CommonEnum::NodeType::ALPHA_CARBON;
}

Network::~Network(){

}

int Network::AddNode(
  int            amino_id,
  string         amino_name,
  vector<int>    atom_ids,
  vector<string> atom_names,
  vector<double> atom_pos_x,
  vector<double> atom_pos_y,
  vector<double> atom_pos_z,
  vector<double> atom_q)
{

    if( !calculateNodeProperties(atom_names, atom_pos_x, atom_pos_y, atom_pos_z, atom_q) ) {
        cout << "Unable to construct node for residual #" << amino_id << " (" << amino_name << ")" << endl;
        return -1;
    }

    mNodes.push_back(
      Node(
        static_cast<int>(mNodes.size()),
        amino_id,
        amino_name,
        mPosX,
        mPosY,
        mPosZ,
        mEstQ,
        atom_ids,
        atom_names,
        atom_pos_x,
        atom_pos_y,
        atom_pos_z,
        atom_q));

    return static_cast<int>(mNodes.size() - 1);
}

bool Network::calculateNodeProperties(vector<string>& atom_names, vector<double>& atom_pos_x, vector<double>& atom_pos_y, vector<double>& atom_pos_z, vector<double>& atom_q) {

  // Common block initialization for all methods
  mWeight = 0.0;
  mPosX   = 0.0;
  mPosY   = 0.0;
  mPosZ   = 0.0;
  mEstQ   = 0.0;

  switch( mEstimateNodeType ) {
    case CommonEnum::NodeType::ALPHA_CARBON :
    default :
      for(int cnt = 0; cnt < static_cast<int>(atom_names.size()); ++cnt) {
        if(atom_names[cnt].compare("CA") == 0) {
          mWeight += 1.0;
          mPosX  += atom_pos_x[cnt];
          mPosY  += atom_pos_y[cnt];
          mPosZ  += atom_pos_z[cnt];
          mEstQ  += atom_q[cnt];
        }
      }
  }

  // Common block to average results for all methods
  if(mWeight == 0.0)
    return false;

  mPosX /= mWeight;
  mPosY /= mWeight;
  mPosZ /= mWeight;
  mEstQ /= mWeight;

  return true;
}

void Network::CreateBackboneConnection(int index1, int index2) {
  if(index1 < 0 || index1 >= static_cast<int>(mNodes.size()) || index2 < 0 || index2 >= static_cast<int>(mNodes.size())) {
    cout << "Network::CreateBackboneConnection - Invalid index" << endl;
    return;
  }

  addConnection(CommonEnum::ConnectionType::SPRING_LEVEL_1, mNodes[index1], mNodes[index2]);
}


void Network::addConnection(CommonEnum::ConnectionType connType, Node& n1, Node& n2) {
  map<CommonEnum::ConnectionType, vector<Connection>>::iterator it = mConnections.find(connType);

  if(it == mConnections.end()) {
    vector<Connection> tmpConn;
    tmpConn.push_back( Connection(n1, n2, connType ) );
    mConnections.insert( pair<CommonEnum::ConnectionType, vector<Connection>>(connType, tmpConn) );
    return;
  }

  it->second.push_back( Connection(n1, n2, connType) );
}


void Network::IdentifyContacts( double cutoff ) {

}


void Network::Print() {
  for(auto node : mNodes)
    node.Print();
/*
  for(map<CommonEnum::ConnectionType, vector<Connection>>::iterator it = mConnections.begin(); it != mConnections.end(); ++it)
    for(auto conn : it->second)
      conn.Print();
      */
}
