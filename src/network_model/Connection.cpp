#include <map>
#include <limits>
#include <iostream>
#include <cmath>

#include "network_model/Connection.h"

using namespace std;



Connection::Connection(Node& n1, Node& n2, CommonEnum::ConnectionType ct) :
mNode1(n1), mNode2(n2), mConnectionType(ct)
{
  mSeparation      = n1.GetSeparationVector(n2);
  mSquaredDistance = mSeparation[0]*mSeparation[0]
                   + mSeparation[1]*mSeparation[1]
                   + mSeparation[2]*mSeparation[2];

  mNode1.AddConnection(ct, n2.GetID());
  mNode2.AddConnection(ct, n1.GetID());
}

Connection::~Connection(){
}


void Connection::Print() {

  cout << "Node IDs: " << mNode1.GetID() << ", " << mNode2.GetID() << endl;
  cout << "\tConnection type: " << mConnectionType << endl;
  cout << "\tSeparation     : ( " << mSeparation[0] << ", " << mSeparation[1] << ", " << mSeparation[2] << " )" << endl;
  cout << "\tSquare Distance: " << mSquaredDistance << endl;
}
