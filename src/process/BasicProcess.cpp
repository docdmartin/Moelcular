#include <iostream>
#include <fstream>

#include <dirent.h>

#include "process/BasicProcess.h"
#include "util/CSVRead.h"

using namespace std;

BasicProcess::BasicProcess(string filename, string pathname) :
  mConfigurationFile( filename ),
  mPathName(pathname),
  mNetworkModel(mParameter)
{
    mMaxIndex = 0;

    // Step #1: Configuration
    loadConfigurationFile();

    // Step #2: Construction
    loadNetworkModel();

    // Step #3: Solve
    mNetworkModel.ConstructLinearResponse();
    //   Build E vector and supply it to the linear solver as the constant vector (b)
    //   Supply omega value and solve
    //     Contract solution with E vector to get get potential
    //     Build T tensor and solve remainder of problem

    mNetworkModel.Print();
}


BasicProcess::~BasicProcess(){

}


void BasicProcess::loadConfigurationFile(){

    CSVRead config_file;

    if(!config_file.OpenCSVFile(mConfigurationFile))
        throw "Unable to open configuration file";

	  vector<string> config_columns;
	  int param_name_index  = static_cast<int>(config_columns.size()); config_columns.push_back("ParameterName" );
	  int param_value_index = static_cast<int>(config_columns.size()); config_columns.push_back("ParameterValue");

	  if(!config_file.ReadCSVHeader(config_columns)) {
		    throw "Configuration file didn't have correct column headers";
	  }

	  while(true) {
		    vector<string> fields = config_file.ReadCSVRecord();
		    if(fields.size() == 0) {
			      break;
		    }

        if(fields[param_name_index].compare("ProteinFolderName") == 0)
            mInputFolder = mPathName + fields[param_value_index];
        if(fields[param_name_index].compare("MaxSpringLength") == 0) {
            double maxSpringLength = atof( fields[param_value_index].c_str() );
            mNetworkModel.SetMaxSpringLength( maxSpringLength );
        }
        if(fields[param_name_index].compare("NodeType") == 0) {
            if(fields[param_value_index].compare("alpha_carbon") == 0) {
                mParameter.SetNodeType( CommonType::NodeType::ALPHA_CARBON );
            }
            else if(fields[param_value_index].compare("mass_weighted") == 0) {
                mParameter.SetNodeType( CommonType::NodeType::MASS_WEIGHTED_MEAN );
            }
        }

        cout << "Parameter: " << fields[param_name_index] << ", has value: " << fields[param_value_index] << endl;
	  }
}


void BasicProcess::loadNetworkModel() {

  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir (mInputFolder.c_str())) != NULL) {

    size_t botDirPos = mInputFolder.find_last_of("/");
    string proteinName = mInputFolder.substr(botDirPos+1, mInputFolder.size() - botDirPos);

    cout << "Protein name: " << proteinName << endl;

    /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      string tmpFileName = ent->d_name;

      if (tmpFileName.compare(0, proteinName.size(), proteinName))
        continue;

      string chain = tmpFileName.substr(proteinName.size()+1, tmpFileName.size() - proteinName.size() - 5);

      countAminoAcid(mInputFolder + "/" + tmpFileName, chain);
    }
    closedir (dir);

    /* print all the files and directories within directory */
    dir = opendir (mInputFolder.c_str());
    while ((ent = readdir (dir)) != NULL) {
      string tmpFileName = ent->d_name;

      if (tmpFileName.compare(0, proteinName.size(), proteinName))
        continue;

      string chain = tmpFileName.substr(proteinName.size()+1, tmpFileName.size() - proteinName.size() - 5);

      loadAminoAcid(mInputFolder + "/" + tmpFileName, chain);
    }
    closedir (dir);

  }

  mNetworkModel.IdentifyContacts();
}

void BasicProcess::countAminoAcid(string input_file, string chain_name) {

  cout << "Loading " << chain_name << " from file: " << input_file << endl;

  CSVRead network_file;

  if(!network_file.OpenCSVFile(input_file))
  throw "Unable to open network file";

  vector<string> network_columns;
  int name_index   = static_cast<int>(network_columns.size()); network_columns.push_back("ATOMNAME" );

  if(!network_file.ReadCSVHeader(network_columns)) {
    throw "Network file didn't have correct column headers";
  }

  while(true) {
    vector<string> fields = network_file.ReadCSVRecord();
    if(fields.size() == 0) {
      break;
    }

    // If atom is not alpha carbon then go to next line
    if(fields[name_index].compare("CA") != 0)
      continue;

    ++mMaxIndex;
  }

}

void BasicProcess::loadAminoAcid(string input_file, string chain_name) {

  cout << "Loading " << chain_name << " from file: " << input_file << endl;

  CSVRead network_file;

  if(mMaxIndex <= 0)
  throw "No nodes";

  cout << input_file << " contains a total of " << mMaxIndex << " alpha carbons." << endl;
  mNetworkModel.AllocateNodes(mMaxIndex);

  vector<string> network_columns;
  int name_index     = static_cast<int>(network_columns.size()); network_columns.push_back("ATOMNAME" );
  int x_index        = static_cast<int>(network_columns.size()); network_columns.push_back("X");
  int y_index        = static_cast<int>(network_columns.size()); network_columns.push_back("Y");
  int z_index        = static_cast<int>(network_columns.size()); network_columns.push_back("Z");
  int q_index        = static_cast<int>(network_columns.size()); network_columns.push_back("Q");
  int atom_id_index  = static_cast<int>(network_columns.size()); network_columns.push_back("ID");
  int resid_index    = static_cast<int>(network_columns.size()); network_columns.push_back("RESID");
  int res_name_index = static_cast<int>(network_columns.size()); network_columns.push_back("RESNAME");

  if(!network_file.OpenCSVFile(input_file))
  throw "Unable to open molecular file";

  if(!network_file.ReadCSVHeader(network_columns)) {
    throw "Molecular file didn't have correct column headers";
  }

  int            resid  = -1;
  string         resname;
  vector<int>    atomid;
  vector<string> atomname;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> q;

  atomid.clear();

  int node_index = -1;
  while(true) {
    vector<string> fields = network_file.ReadCSVRecord();
    if(fields.size() == 0) {
      break;
    }

    if( atoi(fields[resid_index].c_str()) == resid ) {
      atomid  .push_back( atoi(fields[atom_id_index].c_str()) );
      atomname.push_back(      fields[name_index   ]          );
      x       .push_back( atof(fields[x_index      ].c_str()) );
      y       .push_back( atof(fields[y_index      ].c_str()) );
      z       .push_back( atof(fields[z_index      ].c_str()) );
      q       .push_back( atof(fields[q_index      ].c_str()) );

      continue;
    }

    if(atomid.size() > 0) {
      int prev_node_index = node_index;
      node_index          = mNetworkModel.AddNode(resid, resname, atomid, atomname, x, y, z, q);

      if(node_index < 0)
        throw "Process terminated due to error";

      if(prev_node_index >= 0)
        mNetworkModel.CreateBackboneConnection(prev_node_index, node_index);
    }

    resid   = atoi(fields[resid_index].c_str());
    resname = fields[res_name_index];

    atomid  .clear(); atomid  .push_back( atoi(fields[atom_id_index].c_str()) );
    atomname.clear(); atomname.push_back(      fields[name_index   ]          );
    x       .clear(); x       .push_back( atof(fields[x_index      ].c_str()) );
    y       .clear(); y       .push_back( atof(fields[y_index      ].c_str()) );
    z       .clear(); z       .push_back( atof(fields[z_index      ].c_str()) );
    q       .clear(); q       .push_back( atof(fields[q_index      ].c_str()) );
  }

  // Convert last amino acid into a node
  if(atomid.size() > 0) {
    int prev_node_index = node_index;
    node_index          = mNetworkModel.AddNode(resid, resname, atomid, atomname, x, y, z, q);

    if(prev_node_index >= 0)
      mNetworkModel.CreateBackboneConnection(prev_node_index, node_index);
  }

}
