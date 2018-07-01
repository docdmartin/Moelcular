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

}

Network::~Network(){

}

void Network::AddNode(double x, double y, double z) {
    mNodes.push_back(Node(static_cast<int>(mNodes.size()), x, y, z));

    if(mNodes.size() > 1)
      AddConnection(CommonEnum::ConnectionType::SPRING_LEVEL_1, mNodes[mNodes.size()-1], mNodes[mNodes.size()-2]);

    if(mNodes.size() > 2)
      AddConnection(CommonEnum::ConnectionType::SPRING_LEVEL_2, mNodes[mNodes.size()-1], mNodes[mNodes.size()-3]);
}


void Network::AddConnection(CommonEnum::ConnectionType connType, Node& n1, Node& n2) {
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

  for(map<CommonEnum::ConnectionType, vector<Connection>>::iterator it = mConnections.begin(); it != mConnections.end(); ++it)
    for(auto conn : it->second)
      conn.Print();
}
