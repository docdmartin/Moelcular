#include <map>
#include <limits>
#include <iostream>
#include <cmath>

#include "network_model/Connection.h"
#include "util/Common.h"

using namespace std;



Connection::Connection(Common& common_param, Node& n1, Node& n2, CommonType::ConnectionType ct) :
mNode1(n1), mNode2(n2), mConnectionType(ct), mParameter(common_param)
{
  mSpringConstant  = mParameter.mConnectionStrength[mConnectionType];

  // mSeparation is vector pointing from mNode1 to mNode2
  mSeparation      = mNode1.GetSeparationVector(mNode2);

  // mSquaredDistance length of mSeparation squared
  mSquaredDistance = mSeparation[0]*mSeparation[0]
                   + mSeparation[1]*mSeparation[1]
                   + mSeparation[2]*mSeparation[2];

  // Notify nodes of new connection
  mNode1.AddConnection(ct, mNode2.GetID());
  mNode2.AddConnection(ct, mNode1.GetID());
}

Connection::~Connection(){
}


void Connection::Print() {

  cout << "Node IDs: " << mNode1.GetID() << ", " << mNode2.GetID() << endl;
  cout << "\tConnection type: " << mConnectionType << endl;
  cout << "\tSeparation     : ( " << mSeparation[0] << ", " << mSeparation[1] << ", " << mSeparation[2] << " )" << endl;
  cout << "\tSquare Distance: " << mSquaredDistance << endl;
}
