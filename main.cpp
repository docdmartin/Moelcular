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
#include "CSVRead.h"
#include "Solver.h"

using namespace std;


int main(int argc,char **argv)
{
  (void) argc;
  (void) argv;
  try {

    vector<Node> node_array;
    vector<Atom> atom_array;
    map< pair<size_t, size_t>, double > link_map;

    string input_file = "./lys_P1.csv";
    {
      CSVRead network_file;

      if(!network_file.OpenCSVFile(input_file))
        throw "Unable to open network file";

      vector<string> network_columns;
      int atom_name_index   = static_cast<int>(network_columns.size()); network_columns.push_back("ATOMNAME" );
      //int node_name_index   = static_cast<int>(network_columns.size()); network_columns.push_back("RESNAME"  );
      if(!network_file.ReadCSVHeader(network_columns)) {
          throw "Network file didn't have correct column headers";
      }

      int atom_cnt = 0;
      int node_cnt = 0;
      while(true) {
          vector<string> fields = network_file.ReadCSVRecord();
          if(fields.size() == 0) {
              break;
          }

          // If atom is not alpha carbon then go to next line
          if(fields[atom_name_index].compare("CA") == 0)
              ++node_cnt;

          ++atom_cnt;
      }

      node_array.reserve(node_cnt);
      atom_array.reserve(atom_cnt);
    }

    {
      CSVRead network_file;

      if(!network_file.OpenCSVFile(input_file))
          throw "Unable to open network file";

      vector<string> network_columns;
      int atom_name_index   = static_cast<int>(network_columns.size()); network_columns.push_back("ATOMNAME" );
      int node_name_index   = static_cast<int>(network_columns.size()); network_columns.push_back("RESNAME"  );
      int node_id_index     = static_cast<int>(network_columns.size()); network_columns.push_back("RESID"  );
      int x_index           = static_cast<int>(network_columns.size()); network_columns.push_back("X"  );
      int y_index           = static_cast<int>(network_columns.size()); network_columns.push_back("Y"  );
      int z_index           = static_cast<int>(network_columns.size()); network_columns.push_back("Z"  );
      int q_index           = static_cast<int>(network_columns.size()); network_columns.push_back("Q"  );
      if(!network_file.ReadCSVHeader(network_columns)) {
          throw "Network file didn't have correct column headers";
      }

      while(true) {
          vector<string> fields = network_file.ReadCSVRecord();
          if(fields.size() == 0) {
              break;
          }

          atom_array.push_back( Atom() );
          atom_array.back().SetPosition(
            atof( fields[x_index].c_str() ),
            atof( fields[y_index].c_str() ),
            atof( fields[z_index].c_str() ) );
          atom_array.back().SetCharge(
            atof( fields[q_index].c_str() ) );
          atom_array.back().SetName( fields[atom_name_index] );

          if( node_array.empty() || fields[node_id_index] != node_array.back().GetId() ){
            node_array.push_back( Node() );
            node_array.back().SetName( fields[node_name_index] );
            node_array.back().SetId  ( fields[node_id_index] );
          }
          if( node_array.size() > 1 ){
            pair<size_t, size_t> link(node_array.size()-2, node_array.size()-1);
            link_map.insert( pair<pair<size_t, size_t>, double>(link, 100.0) );

            // quick check for near nodes, link strength 1.0
            // later this will be done in a hash
            double* node_pos = node_array.back().GetPosition();
            double* tmp_node;
            double dx, dy, dz;
            for(size_t node_index = 0; node_index < node_array.size()-2; ++node_index){
              tmp_node = node_array[node_index].GetPosition();
              dx = *(node_pos  ) - *(tmp_node  );
              dy = *(node_pos+1) - *(tmp_node+1);
              dz = *(node_pos+2) - *(tmp_node+2);
              if( fabs(dx) > 15.0 || fabs(dy) > 15.0 || fabs(dz) > 15.0 )
                continue;

              if( dx*dx + dy*dy + dz*dz < 225. ){
                pair<size_t, size_t> link(node_index, node_array.size()-1);
                link_map.insert( pair<pair<size_t, size_t>, double>(link, 1.0) );
              }
            }
          }

          node_array.back().AddAtom( &(atom_array.back()) );

          // If atom is not alpha carbon then go to next line
          if( ElementType::C == atom_array.back().GetName() && atom_array.back().GetGreek() == 'A' )
              node_array.back().SetPosition(
                atof( fields[x_index].c_str() ),
                atof( fields[y_index].c_str() ),
                atof( fields[z_index].c_str() ) );
      }
    }

    vector<Matrix3x3> hessian;
    hessian.reserve(link_map.size());
    map< pair<size_t, size_t>, double >::iterator it = link_map.begin();
    for(; it != link_map.end(); ++it){
      hessian.push_back( Matrix3x3(it->first.first, it->first.second,
        node_array[it->first.first ].GetPosition(),
        node_array[it->first.second].GetPosition(), it->second) );
    }

    // given a probe position generate electric potential
    double  probe_pos[3] = {-14.167, 4.687, -0.717};
    double* electric_potential = new double[3*node_array.size()]();

    size_t node_index = 0;
    for( Node& node : node_array ){
      node.SetElectricPotential( &(electric_potential[node_index]), &(probe_pos[0]) );
      node_index += 3;
    }

    double freq = 1.0;
    vector<double> frequency;
    for(int cnt = 0; cnt < 5; ++cnt){
      frequency.push_back( freq );
      freq /= 2.;
    }

    Solver solver( 3*node_array.size(), hessian, frequency, electric_potential );

    vector<pair<double, double>> response = solver.CalculateResponse();

  //  for( Node& node : node_array ){
  //    node.Print();
  //  }
    cout << "Total number of nodes: " << node_array.size() << endl;
    cout << "Total number of links: " << link_map.size() << endl;

    for( size_t cnt = 0; cnt < frequency.size(); ++cnt ){
      cout << "Frquency " << frequency[cnt] << " has response: " << response[cnt].first << " + " << response[cnt].second << " i" << endl;
    }



    delete[] electric_potential;
  }
  catch(string e) {
	  cout << "Terminating due to error: " << e << endl;
	  return 1;
  }

  return 0;
}
