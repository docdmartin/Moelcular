#include <iomanip>

#include "network_model/Node.h"

Node::Node(Common& common_param,
  int                  id,
  int                  amino_id,
  string               amino_name,
  double               node_pos_x,
  double               node_pos_y,
  double               node_pos_z,
  double               node_q,
  vector<int>          atom_ids,
  vector<string>       atom_names,
  vector<double>       atom_pos_x,
  vector<double>       atom_pos_y,
  vector<double>       atom_pos_z,
  vector<double>       atom_q)
  : mNodeID(id), mParameter(common_param)
{

    mPos[0] =  node_pos_x;
    mPos[1] =  node_pos_y;
    mPos[2] =  node_pos_z;
    mCharge =  node_q;

    mResidualId   = amino_id;
    mResidualName = amino_name;

    mAtomId     = atom_ids;
    mAtomName   = atom_names;
    mAtomPosX   = atom_pos_x;
    mAtomPosY   = atom_pos_y;
    mAtomPosZ   = atom_pos_z;
    mAtomCharge = atom_q;
}

Node::~Node(){
}

void Node::SetPosition(double x, double y, double z) {
    mPos[0] = x;
    mPos[1] = y;
    mPos[2] = z;
}

vector<double> Node::GetSeparationVector(Node& tmpNode) {
  // Seperation between two vectors in Cartesian coordinates
  vector<double> los(3, 0.0);
  los[0] = tmpNode.GetX() - mPos[0];
  los[1] = tmpNode.GetY() - mPos[1];
  los[2] = tmpNode.GetZ() - mPos[2];

  return los;
}

void Node::AddConnection(CommonType::ConnectionType ct, int id) {
  map< CommonType::ConnectionType, set<int> >::iterator it = mConnections.find(ct);
  if(it == mConnections.end()) {
    set<int> tmpSet;
    tmpSet.insert( id );
    mConnections.insert( pair< CommonType::ConnectionType, set<int> >(ct, tmpSet) );
    return;
  }

  it->second.insert( id );
}

bool Node::IsConnected(int id, CommonType::ConnectionType ct) {

  if(ct == CommonType::ConnectionType::ALL_CONNECTION) {
    for ( int enumInt = CommonType::ConnectionType::NO_CONNECTION + 1; enumInt != CommonType::ConnectionType::ALL_CONNECTION; ++enumInt ) {
      bool tmp = IsConnected(id, static_cast<CommonType::ConnectionType>(enumInt));
      if(tmp)
        return true;
    }
    return false;
  }

  map< CommonType::ConnectionType, set<int> >::iterator it = mConnections.find(ct);
  if(it == mConnections.end())
    return false;

  if(it->second.find(id) == it->second.end())
    return false;

  return true;
}


void Node::Print() {
  cout << "Node ID: " << mNodeID << ", type: " << mParameter.mNodeTypeDef[mParameter.GetNodeType()] << endl;
  cout << "\tPosition : ( " << mPos[0] << ", " << mPos[1] << ", " << mPos[2] << " )" << endl;
  cout << "\tCharge   : " << mCharge << endl;
  cout << "\tResidual : #" << mResidualId << " ( " << mResidualName << " )" << endl;
  cout << "\tAtoms    : #" << setw(4) << mAtomId[0]
                           << " ( " << setw(4) << mAtomName[0]
                           << " ), Position: ( " << setw(6) << setprecision(4) << mAtomPosX[0]
                           << ", " << setw(6) << setprecision(4) << mAtomPosY[0]
                           << ", " << setw(6) << setprecision(4) << mAtomPosZ[0]
                           << " ) and Charge : " << setw(6) << setprecision(4) << mAtomCharge[0] << endl;
  for(int cnt = 1; cnt < static_cast<int>(mAtomId.size()); ++cnt)
    cout << "\t           #"  << setw(4) << mAtomId[cnt]
                              << " ( "  << setw(4) << mAtomName[cnt]
                              << " ), Position: ( " << setw(6) << setprecision(4) << mAtomPosX[cnt]
                              << ", " << setw(6) << setprecision(4) << mAtomPosY[cnt]
                              << ", " << setw(6) << setprecision(4) << mAtomPosZ[cnt]
                              << " ) and Charge : " << setw(6) << setprecision(4) << mAtomCharge[cnt] << endl;

  for(map< CommonType::ConnectionType, set<int> >::iterator it = mConnections.begin(); it != mConnections.end(); ++it) {
    cout << "\tConnection type: " << it->first << endl;
    cout << "\t";
    for(int tmp : it->second)
      cout << "\t" << tmp;
    cout << endl;
  }
}
