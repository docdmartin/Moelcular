#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include "dirent.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Common.h"
#include "GlobalMaps.h"
#include "EnumToString.h"
#include "StringToEnum.h"


#include "Atom.h"
#include "Node.h"
#include "Matrix3x3.h"
#include "Kernel.h"
#include "CSVRead.h"
#include "Solver.h"
#include "SQLite.h"
#include "Manager.h"

using namespace std;


int main(int argc,char **argv)
{
  (void) argc;
  (void) argv;
  try {

    Manager cnfg_manager;
    string data_base = "../Test/create_SQLITE_DB/lys_P1.db";
    SQLite sql_lite( data_base, cnfg_manager );

    cnfg_manager.BuildNodes();

    vector<double>                     frequency          = cnfg_manager.GetFrequency();
    vector<Node>&                      node_array         = cnfg_manager.GetNodeArray();
    map<pair<size_t, size_t>, double>& link_map           = cnfg_manager.GetLinkMap();
    map< int, vector<double> >&        probe_location_map = cnfg_manager.GetProbe();
    vector< pair<double, double> >&    friction_kernel    = cnfg_manager.GetKernel();

    Kernel kernel( node_array );

    vector<Matrix3x3> hessian;
    hessian.reserve(link_map.size());
    map< pair<size_t, size_t>, double >::iterator it = link_map.begin();
    for(; it != link_map.end(); ++it){
      hessian.push_back( Matrix3x3(it->first.first, it->first.second,
        node_array[it->first.first ].GetPosition(),
        node_array[it->first.second].GetPosition(), it->second) );
    }

    for( map< int, vector<double> >::iterator it = probe_location_map.begin(); it != probe_location_map.end(); ++it ){

      cout << "Solving for probe id: " << it->first << ", ( " << it->second[0] << ", " << it->second[1] << ", " << it->second[2] << " )" << endl;

      double  probe_pos[3] = {it->second[0], it->second[1], it->second[2]};
      double* electric_potential = new double[3*node_array.size()]();

      size_t node_index = 0;
      for( Node& node : node_array ){
        node.SetElectricPotential( &(electric_potential[node_index]), &(probe_pos[0]) );
        node_index += 3;
      }

      kernel.RemoveProjection( electric_potential );

      Solver solver( 3*node_array.size(), hessian, frequency, electric_potential, kernel );//

      double total_weight = 0.;
      vector<pair<double, double>> response;
      for( auto& fk : friction_kernel ){

        solver.Reset();

        cout << "Kernel: " << fk.first << ", " << fk.second << endl;
        total_weight += fk.first;

        vector<pair<double, double>> rspn = solver.CalculateResponse( fk.second );
        if( response.empty() ){
          response = rspn;
          for( pair<double, double>& r : response ){
            r.first  *= fk.first;
            r.second *= fk.first;
          }
        }
        else{
          for( size_t cnt = 0; cnt < rspn.size(); ++cnt ){
            response[cnt].first  += rspn[cnt].first  * fk.first;
            response[cnt].second += rspn[cnt].second * fk.first;
          }
        }
      }
      for( size_t cnt = 0; cnt < response.size(); ++cnt ){
        response[cnt].first  /= total_weight;
        response[cnt].second /= total_weight;
      }

      SQLite sql_write( data_base, it->first, frequency, response );

      delete[] electric_potential;
    }
  }
  catch(string e) {
	  cout << "Terminating due to error: " << e << endl;
	  return 1;
  }

  return 0;
}
