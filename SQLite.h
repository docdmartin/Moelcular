#ifndef SQLITE_H
#define SQLITE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include "Manager.h"

using namespace std;

class SQLite {

public:

  SQLite( string& database_name, int probe_id, vector<double>& freq, vector< pair<double, double> >& response );

  SQLite( string& database_name, Manager& mngr );
  ~SQLite(){}

private:

  static int sql_config    ( void* data, int argc, char** argv, char** col_name );
  static int sql_count_node( void* data, int argc, char** argv, char** col_name );
  static int sql_count_atom( void* data, int argc, char** argv, char** col_name );
  static int sql_get_atoms ( void* data, int argc, char** argv, char** col_name );
  static int sql_get_probe ( void* data, int argc, char** argv, char** col_name );
  static int sql_get_kernel( void* data, int argc, char** argv, char** col_name );

};

#endif
