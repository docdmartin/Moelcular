#include <iostream>
#include <fstream>



#include "process/BasicProcess.h"
#include "util/CSVRead.h"

using namespace std;

BasicProcess::BasicProcess(string filename) :
  mConfigurationFile( filename )
{

    loadConfigurationFile();

    loadNetworkModel();

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

        if(fields[param_name_index].compare("MolecularFileName") == 0)
            mInputFile = fields[param_value_index];

        cout << "Parameter: " << fields[param_name_index] << ", has value: " << fields[param_value_index] << endl;
	}
}


void BasicProcess::loadNetworkModel() {

  CSVRead network_file;

  /*=================================================
  */
  if(!network_file.OpenCSVFile(mInputFile))
  throw "Unable to open network file";

  vector<string> network_columns;
  int x_index   = static_cast<int>(network_columns.size()); network_columns.push_back("X");

  if(!network_file.ReadCSVHeader(network_columns)) {
    throw "Network file didn't have correct column headers";
  }

  int max_index = 0;
  while(true) {
    vector<string> fields = network_file.ReadCSVRecord();
    if(fields.size() == 0) {
      break;
    }

    ++max_index;
  }

  if(max_index <= 0)
  throw "No nodes";

  cout << mInputFile << " contains a total of " << max_index << " nodes." << endl;
  mNetworkModel.AllocateNodes(max_index);
  /*==================================================
  */


  /*=================================================
  */
  if(!network_file.OpenCSVFile(mInputFile)) /* DAN - WHY IS THIS HERE? ISN'T THIS DONE ABOVE ALREADY ONCE? (2018-07-14)*/
  throw "Unable to open molecular file";

  int y_index       = static_cast<int>(network_columns.size()); network_columns.push_back("Y");
  int z_index       = static_cast<int>(network_columns.size()); network_columns.push_back("Z");

  if(!network_file.ReadCSVHeader(network_columns)) {
    throw "Molecular file didn't have correct column headers";
  }

  while(true) {
    vector<string> fields = network_file.ReadCSVRecord();
    if(fields.size() == 0) {
      break;
    }

    mNetworkModel.AddNode(
      atof(fields[x_index].c_str()),
      atof(fields[y_index].c_str()),
      atof(fields[z_index].c_str()) );
    }

    mNetworkModel.Print();
  }
