#include <map>
#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "ElasticNetworkModel.h"

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



ElasticNetworkModel::ElasticNetworkModel() {
    mMaxSpringLength = 1.0;
    mMaxRangeSq      = mMaxSpringLength * mMaxSpringLength;
}

ElasticNetworkModel::~ElasticNetworkModel(){

}

int ElasticNetworkModel::AddReferencePoint( double x, double y, double z ) {
    int index = static_cast<int>( mReferencePoints.size() );
    mReferencePoints.push_back( ReferencePoint(x, y, z, mNodes, mHessianMatrix) );

    return index;
}

vector< pair<double, double> > ElasticNetworkModel::SingleModeFrequencyResponse(int ref_index, double omega) {
    vector< pair<double, double> > complex_potential;
    if( ref_index < 0 || ref_index >= static_cast<int>( mReferencePoints.size() ) )
        return complex_potential;

    complex_potential = mReferencePoints[ref_index].SingleModeFrequencyResponse(omega);
    return complex_potential;
}
vector< pair<double, double> > ElasticNetworkModel::DualModeFrequencyResponse(int ref_index, double omega1, double omega2, double a) {
    vector< pair<double, double> > complex_potential;
    if( ref_index < 0 || ref_index >= static_cast<int>( mReferencePoints.size() ) )
        return complex_potential;

    complex_potential = mReferencePoints[ref_index].DualModeFrequencyResponse(omega1, omega2, a);
    return complex_potential;
}

void ElasticNetworkModel::SetNodeType(string node_type) {
    if(node_type.compare("alpha_carbon") == 0) {
        mParameter.SetNodeType( CommonType::NodeType::ALPHA_CARBON );
    }
    else if(node_type.compare("mass_weighted") == 0) {
        mParameter.SetNodeType( CommonType::NodeType::MASS_WEIGHTED_MEAN );
    }
}

int ElasticNetworkModel::AddNode(
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

    voxelizeNode();

    mNodes.push_back(
        Node(
            mParameter,
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
            atom_q
        )
    );

    return static_cast<int>(mNodes.size() - 1);
}

bool ElasticNetworkModel::calculateNodeProperties(
    vector<string>& atom_names,
    vector<double>& atom_pos_x,
    vector<double>& atom_pos_y,
    vector<double>& atom_pos_z,
    vector<double>& atom_q)
{

    // Common block initialization for all methods
    mWeight = 0.0;
    mPosX   = 0.0;
    mPosY   = 0.0;
    mPosZ   = 0.0;
    mEstQ   = 0.0;

    switch( mParameter.GetNodeType() ) {

      // Used when nodes are to be constructed according to alpha carbon
      case CommonType::NodeType::ALPHA_CARBON :

        for(int cnt = 0; cnt < static_cast<int>(atom_names.size()); ++cnt) {
            if(atom_names[cnt].compare("CA") == 0) {
                mWeight += 1.0;
                mPosX  += atom_pos_x[cnt];
                mPosY  += atom_pos_y[cnt];
                mPosZ  += atom_pos_z[cnt];
            }
            mEstQ  += atom_q[cnt];
        }
        break;


      // Used when nodes are to be constructed according to the weighted mean
      case CommonType::NodeType::MASS_WEIGHTED_MEAN :

        for(int cnt = 0; cnt < static_cast<int>(atom_names.size()); ++cnt) {
            string abbrev = atom_names[cnt].substr(0,1);
            map<string, CommonType::ElementType>::iterator it = mKnownElements.find(abbrev);
            if(it == mKnownElements.end()) {
                map<CommonType::ElementType, Common::ElementData>::iterator iter = mParameter.mPeriodicTable.begin();
                for(; iter != mParameter.mPeriodicTable.end(); ++iter) {
                    if(abbrev.compare(iter->second.abbreviation) == 0) {
                        mKnownElements.insert( pair<string, CommonType::ElementType>(iter->second.abbreviation, iter->first) );
                        it = mKnownElements.find(abbrev);
                        break;
                    }
                }
            }

            if(it == mKnownElements.end()) {
                cout << "Unable to find element: " << atom_names[cnt] << endl;
                continue;
            }

            mWeight += mParameter.mPeriodicTable[it->second].mass;
            mPosX   += mParameter.mPeriodicTable[it->second].mass * atom_pos_x[cnt];
            mPosY   += mParameter.mPeriodicTable[it->second].mass * atom_pos_y[cnt];
            mPosZ   += mParameter.mPeriodicTable[it->second].mass * atom_pos_z[cnt];
            mEstQ   += atom_q[cnt];
        }

        break;


      // Configuration file failed to identify how nodes are to be constructed
      // Please insert one of the following lines into configuration file
      //  NodeType,alpha_carbon
      //  NodeType,mass_weighted
      case CommonType::NodeType::UNDEFINED :
      default :
        cout << "Node type " << mParameter.GetNodeType() << " is undefined" << endl;
        mWeight = 0.0;
        break;
    }

    // Remove numerical noise built into the charge
    mEstQ = 1e-4 * floor(mEstQ * 1e4 + 0.5);


    // Common block to average results for all methods
    if(mWeight == 0.0)
        return false;

    mPosX /= mWeight;
    mPosY /= mWeight;
    mPosZ /= mWeight;

    return true;
}

void ElasticNetworkModel::CreateBackboneConnection(int index1, int index2) {
  if(index1 < 0 || index1 >= static_cast<int>(mNodes.size()) || index2 < 0 || index2 >= static_cast<int>(mNodes.size())) {
    cout << "ElasticNetworkModel::CreateBackboneConnection - Invalid index" << endl;
    return;
  }

  addConnection(CommonType::ConnectionType::SPRING_LEVEL_1, mNodes[index1], mNodes[index2]);
}


void ElasticNetworkModel::addConnection(CommonType::ConnectionType connType, Node& n1, Node& n2) {
  map<CommonType::ConnectionType, vector<Connection>>::iterator it = mConnections.find(connType);

  if(it == mConnections.end()) {
    vector<Connection> tmpConn;
    tmpConn.push_back( Connection(mParameter, n1, n2, connType ) );
    mConnections.insert( pair<CommonType::ConnectionType, vector<Connection>>(connType, tmpConn) );
    return;
  }

  it->second.push_back( Connection(mParameter, n1, n2, connType) );
}


void ElasticNetworkModel::IdentifyContacts() {
  map<int, map<int, map<int, vector<int> > > >::iterator xit;
  map<         int, map<int, vector<int> >   >::iterator yit;
  map<                  int, vector<int>     >::iterator zit;

  map<int, map<int, map<int, vector<int> > > >::iterator xiter;
  map<         int, map<int, vector<int> >   >::iterator yiter;
  map<                  int, vector<int>     >::iterator ziter;

  int pri_vector_index;
  for(    xit = mNodeVoxels.begin(); xit != mNodeVoxels.end(); ++xit) {
    for(  yit = xit->second.begin(); yit != xit->second.end(); ++yit) {
      for(zit = yit->second.begin(); zit != yit->second.end(); ++zit) {
        for(int ii = 0; ii < static_cast<int>(zit->second.size()); ++ii) {
          pri_vector_index = zit->second.at(ii);

          // Look in same voxel
          //   only look forward since node behind pri_vector_index have already been tested
          for(int iii = ii+1; iii < static_cast<int>(zit->second.size()); ++iii) {
            testConnection(pri_vector_index, zit->second.at(iii));
          }

          // need to check for same x, y but z + 1
          ziter = yit->second.find( zit->first + 1 );
          if(ziter != yit->second.end()) {
            for(int iii = 0; iii < static_cast<int>(ziter->second.size()); ++iii) {
              testConnection(pri_vector_index, ziter->second.at(iii));
            }
          }

          // need to check for same x but y + 1 and z-1, z, z+1
          yiter = xit->second.find( yit->first + 1 );
          if(yiter != xit->second.end()) {
            for(int z_cnt = -1; z_cnt <= 1; ++z_cnt) {
              ziter = yiter->second.find( zit->first + z_cnt );
              if(ziter == yiter->second.end())
                continue;

              for(int iii = 0; iii < static_cast<int>(ziter->second.size()); ++iii) {
                testConnection(pri_vector_index, ziter->second.at(iii));
              }
            }
          }

          // need to check for x+1, y-1, y, y+1, z-1, z, z+1
          xiter = mNodeVoxels.find( xit->first + 1 );
          if(xiter != mNodeVoxels.end()) {
            for(int y_cnt = -1; y_cnt <= 1; ++y_cnt) {
              yiter = xiter->second.find( yit->first + y_cnt );
              if(yiter == xiter->second.end())
                continue;

              for(int z_cnt = -1; z_cnt <= 1; ++z_cnt) {
                ziter = yiter->second.find( zit->first + z_cnt );
                if(ziter == yiter->second.end())
                  continue;

                for(int iii = 0; iii < static_cast<int>(ziter->second.size()); ++iii) {
                  testConnection(pri_vector_index, ziter->second.at(iii));
                }
              }
            }
          }
        }
      }
    }
  }
}

/*
  All connections are established
  Begin construction of linear response function
*/
void ElasticNetworkModel::ConstructLinearResponse() {

  map<CommonType::ConnectionType, vector<Connection>>::iterator it = mConnections.begin();
  for(; it != mConnections.end(); ++it) {
    for(Connection conn : it->second) {
      mHessianMatrix.SetConnection(conn.GetNodeId2(), conn.GetNodeId1(), conn.GetSeparation(), conn.GetSpringConstant());
    }
  }
}

void ElasticNetworkModel::testConnection(int p_index, int s_index) {
  // if nodes are already linked then continue to next node
  if( mNodes[p_index].IsConnected(s_index, CommonType::ConnectionType::ALL_CONNECTION) )
    return;

  // range suare test
  // mSeparation is vector pointing from mNode1 to mNode2
  vector<double> seperation = mNodes[p_index].GetSeparationVector(mNodes[s_index]);

  // mSquaredDistance length of mSeparation squared
  double squared_distance = seperation[0]*seperation[0]
                   + seperation[1]*seperation[1]
                   + seperation[2]*seperation[2];

  // if pass test then create secondary link
  if(squared_distance < mMaxRangeSq) {
    addConnection(CommonType::ConnectionType::SPRING_LEVEL_2, mNodes[p_index], mNodes[s_index]);
  }
}

void ElasticNetworkModel::voxelizeNode() {
  int x_index = floor( mPosX / mMaxSpringLength );
  int y_index = floor( mPosY / mMaxSpringLength );
  int z_index = floor( mPosZ / mMaxSpringLength );

  int node_index = static_cast<int>( mNodes.size() );

  map<int, map<int, map<int, vector<int> > > >::iterator xit = mNodeVoxels.find(x_index);
  if(xit == mNodeVoxels.end()) {
    vector<int>                       n_index;
    map<int, vector<int> >            zmap;
    map<int, map<int, vector<int> > > ymap;

    n_index.push_back( node_index );
    zmap[z_index]        = n_index;
    ymap[y_index]        = zmap;
    mNodeVoxels[x_index] = ymap;

    return;
  }

  map<int, map<int, vector<int> > >::iterator yit = xit->second.find(y_index);
  if(yit == xit->second.end()) {
    vector<int>                       n_index;
    map<int, vector<int> >            zmap;

    n_index.push_back( node_index );
    zmap[z_index]        = n_index;
    xit->second.insert(pair<int, map<int, vector<int> >  >(y_index, zmap));

    return;
  }

  map<int, vector<int> >::iterator zit = yit->second.find(z_index);
  if(zit == yit->second.end()) {
    vector<int> n_index;

    n_index.push_back( node_index );
    yit->second.insert(pair<int, vector<int> >(z_index, n_index));

    return;
  }

  zit->second.push_back( node_index );
}


void ElasticNetworkModel::Print() {
/*
  map<CommonType::ElementType, Common::ElementData>::iterator it    = mParameter.mPeriodicTable.begin();
  map<CommonType::ElementType, Common::ElementData>::iterator endIt = mParameter.mPeriodicTable.end();
  for(; it != endIt; ++it) {
      cout << setiosflags(ios::right);
      cout << setiosflags(ios::fixed);
      cout << "Element #"    << setw(3) << it->first
           << ", abbr.: "    << std::left  << setw(2) << it->second.abbreviation
           << ", name: "     << setw(13)<< it->second.name
           << ", and mass: " << std::right << setw(10) << setprecision(6) << it->second.mass << endl;
  }
  return;
*/

  cout << "Nodes" << endl;
  for(auto node : mNodes)
    node.Print();

    return;

  cout << "Connections" << endl;
  for(map<CommonType::ConnectionType, vector<Connection>>::iterator it = mConnections.begin(); it != mConnections.end(); ++it)
    for(auto conn : it->second)
      conn.Print();

  mHessianMatrix.PrintMatrixA();

return;

  cout << "Voxels" << endl;
  map<int, map<int, map<int, vector<int> > > >::iterator xit = mNodeVoxels.begin();
  map<         int, map<int, vector<int> >   >::iterator yit;
  map<                  int, vector<int>     >::iterator zit;
  for(; xit != mNodeVoxels.end(); ++xit) {
    for(yit = xit->second.begin(); yit != xit->second.end(); ++yit) {
      for(zit = yit->second.begin(); zit != yit->second.end(); ++zit) {
        cout << "Voxel (" << xit->first << ", " << yit->first << ", " << zit->first << ") has " << zit->second.size() << " nodes" << endl;
      }
    }
  }
}
