
#include "SQLite.h"
#include "sqlite3.h"

using namespace std;

SQLite::SQLite( string& database_name, int probe_id, vector<double>& freq, vector< pair<double, double> >& response ){
  sqlite3* db;
  int exit = 0;
  exit = sqlite3_open( database_name.c_str(), &db );
  if( exit ){
    throw string("Error opening database: " + database_name);
  }

  char buf[512];

  sqlite3_exec( db, "BEGIN TRANSACTION", NULL, NULL, NULL );
  for( size_t cnt = 0; cnt < freq.size(); ++cnt ){
    sprintf(buf, "delete from Response where ProbeId = %d", probe_id);
    sqlite3_exec( db, buf, NULL, NULL, NULL );
  }
  sqlite3_exec( db, "END TRANSACTION", NULL, NULL, NULL );


  sqlite3_exec( db, "BEGIN TRANSACTION", NULL, NULL, NULL );
  for( size_t cnt = 0; cnt < freq.size(); ++cnt ){
    sprintf(buf, "INSERT INTO Response VALUES (%d, %g, %g, %g)", probe_id, freq[cnt], response[cnt].first, response[cnt].second);
    sqlite3_exec( db, buf, NULL, NULL, NULL );
  }
  sqlite3_exec( db, "END TRANSACTION", NULL, NULL, NULL );

  sqlite3_close( db );
}

SQLite::SQLite( string& database_name, Manager& cnfg_manager ){

  sqlite3* db;
  int exit = 0;
  exit = sqlite3_open( database_name.c_str(), &db );
  if( exit ){
    throw string("Error opening database: " + database_name);
  }

  string cnfg_str     = "select * from Configuration";
  string node_cnt_str = "select count(*) nodes from ( select count(*) cnt from PDB_PSF group by segid, resid )";
  string atom_cnt_str = "select count(*) atoms from PDB_PSF";
  string protein_str  = "select segid, resid, atomname, ax, ay, az, charge, mass from PDB_PSF order by segid, resid";
  string probe_str    = "select * from Probe";
  string kernel_str   = "select * from Kernel";

  int return_value;

  return_value = sqlite3_exec(db, cnfg_str.c_str(), &SQLite::sql_config, (void*)&cnfg_manager, NULL);
  if( return_value != SQLITE_OK )
    throw string( "SQLite3 failed to complete - configuration query" );

  return_value = sqlite3_exec(db, node_cnt_str.c_str(), &SQLite::sql_count_node, (void*)&cnfg_manager, NULL);
  if( return_value != SQLITE_OK )
    throw string( "SQLite3 failed to complete - node count query" );

  return_value = sqlite3_exec(db, atom_cnt_str.c_str(), &SQLite::sql_count_atom, (void*)&cnfg_manager, NULL);
  if( return_value != SQLITE_OK )
    throw string( "SQLite3 failed to complete - atom count query" );

  return_value = sqlite3_exec(db, protein_str.c_str(), &SQLite::sql_get_atoms, (void*)&cnfg_manager, NULL);
  if( return_value != SQLITE_OK )
    throw string( "SQLite3 failed to complete - atom query" );

  return_value = sqlite3_exec(db, probe_str.c_str(), &SQLite::sql_get_probe, (void*)&cnfg_manager, NULL);
  if( return_value != SQLITE_OK )
    throw string( "SQLite3 failed to complete - probe query" );

  return_value = sqlite3_exec(db, kernel_str.c_str(), &SQLite::sql_get_kernel, (void*)&cnfg_manager, NULL);
  if( return_value != SQLITE_OK )
    throw string( "SQLite3 failed to complete - kernel query" );

  sqlite3_close( db );

}

int SQLite::sql_get_kernel( void* data, int argc, char** argv, char** col_name ){
  int    id;
  double weight;
  double freq_mult;

  for(int cnt = 0; cnt < argc; ++cnt){
    string tmp = col_name[cnt];
    if( tmp == "id" ){
      id = static_cast<int>( atoi(argv[cnt]) );
    }
    else if( tmp == "weight" ){
      weight = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "freq_mult" ){
      freq_mult = static_cast<double>( atof(argv[cnt]) );
    }
  }

  Manager* const mngr = static_cast<Manager*>(data);
  mngr->AddKernel( id, weight, freq_mult );

  return 0;
}

int SQLite::sql_get_probe( void* data, int argc, char** argv, char** col_name ){
  int    id;
  double pos[3] = {0., 0., 0.};

  for(int cnt = 0; cnt < argc; ++cnt){
    string tmp = col_name[cnt];
    if( tmp == "id" ){
      id = static_cast<int>( atoi(argv[cnt]) );
    }
    else if( tmp == "x" ){
      pos[0] = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "y" ){
      pos[1] = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "z" ){
      pos[2] = static_cast<double>( atof(argv[cnt]) );
    }
  }

  Manager* const mngr = static_cast<Manager*>(data);
  mngr->AddProbe( id, pos );

  return 0;
}

int SQLite::sql_get_atoms( void* data, int argc, char** argv, char** col_name ){
  double pos[3] = {0., 0., 0.};
  double mass   = 0.;
  double charge = 0.;
  string segid;
  string resid;
  string atomname;

  for(int cnt = 0; cnt < argc; ++cnt){
    string tmp = col_name[cnt];
    if( tmp == "segid" ){
      segid = argv[cnt];
    }
    else if( tmp == "resid" ){
      resid = argv[cnt];
    }
    else if( tmp == "atomname" ){
      atomname = argv[cnt];
    }
    else if( tmp == "ax" ){
      pos[0] = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "ay" ){
      pos[1] = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "az" ){
      pos[2] = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "mass" ){
      mass = static_cast<double>( atof(argv[cnt]) );
    }
    else if( tmp == "charge" ){
      charge = static_cast<double>( atof(argv[cnt]) );
    }
  }

  Manager* const mngr = static_cast<Manager*>(data);
  mngr->AddAtom( segid, resid, atomname, pos, mass, charge );

  return 0;
}

int SQLite::sql_count_node( void* data, int argc, char** argv, char** col_name ){
  (void) col_name;
  (void) argc;
  Manager* const mngr = static_cast<Manager*>(data);
  mngr->SetNodeCount( static_cast<int>( atoi(argv[0]) ) );

  return 0;
}

int SQLite::sql_count_atom( void* data, int argc, char** argv, char** col_name ){
  (void) col_name;
  (void) argc;

  Manager* const mngr = static_cast<Manager*>(data);
  mngr->SetAtomCount( static_cast<int>( atoi(argv[0]) ) );

  return 0;
}


int SQLite::sql_config( void* data, int argc, char** argv, char** col_name ){

  Manager* const mngr = static_cast<Manager*>(data);

  int cnfg_id = -1;
  for( int cnt = 0; cnt < argc; ++cnt ){
    string tmp = col_name[cnt];
    if( tmp == "id" ){
      cnfg_id = static_cast<int>( atoi(argv[cnt]) );
    }
  }
  if( cnfg_id < 0 )
    throw string( "Invalid configuration ID" );

  mngr->CreateConfiguration( cnfg_id );

  for( int cnt = 0; cnt < argc; ++cnt ){
    string tmp = col_name[cnt];
    if( tmp == "id" )
      continue;

    mngr->SetParameter( cnfg_id, col_name[cnt], argv[cnt] );
  }

  return 0;
}
