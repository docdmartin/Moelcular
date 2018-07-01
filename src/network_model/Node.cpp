#include "network_model/Node.h"

Node::Node(int id, double x, double y, double z)
  : mNodeID(id)
{
    mPos[0] =  x;
    mPos[1] =  y;
    mPos[2] =  z;
}

Node::~Node(){
}

void Node::SetPosition(double x, double y, double z) {
    mPos[0] = x;
    mPos[1] = y;
    mPos[2] = z;
}

vector<double> Node::GetSeparationVector(Node& tmpNode) {
  vector<double> los;
  los.push_back(tmpNode.GetX() - mPos[0]);
  los.push_back(tmpNode.GetY() - mPos[1]);
  los.push_back(tmpNode.GetZ() - mPos[2]);

  return los;
}

void Node::AddConnection(CommonEnum::ConnectionType ct, int id) {
  map< CommonEnum::ConnectionType, set<int> >::iterator it = mConnections.find(ct);
  if(it == mConnections.end()) {
    set<int> tmpSet;
    tmpSet.insert( id );
    mConnections.insert( pair< CommonEnum::ConnectionType, set<int> >(ct, tmpSet) );
    return;
  }

  it->second.insert( id );
}

bool Node::IsConnected(int id, CommonEnum::ConnectionType ct) {
  map< CommonEnum::ConnectionType, set<int> >::iterator it = mConnections.find(ct);
  if(it == mConnections.end())
    return false;

  if(it->second.find(id) == it->second.end())
    return false;

  return true;
}


void Node::Print() {
  cout << "Node ID: " << mNodeID << ", type: " << mNodeType << endl;
  cout << "\tPosition : ( " << mPos[0] << ", " << mPos[1] << ", " << mPos[2] << " )" << endl;
  for(map< CommonEnum::ConnectionType, set<int> >::iterator it = mConnections.begin(); it != mConnections.end(); ++it) {
    cout << "\tConnection type: " << it->first << endl;
    cout << "\t";
    for(int tmp : it->second)
      cout << "\t" << tmp;
    cout << endl;
  }
}
