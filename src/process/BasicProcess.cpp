#include <iostream>
#include <fstream>



#include "BasicProcess.h"
#include "../../src/util/CSVRead.h"

using namespace std;

BasicProcess::BasicProcess(string filename) :
  mConfigurationFile( filename ),
  mSpringCutoffLength( 1.0 )
{
    
    loadConfigurationFile();
    
    loadMolecule();
    
    mMolecule.IdentifyContacts( mSpringCutoffLength );
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
            mMoleculeFile = fields[param_value_index];
        else if(fields[param_name_index].compare("MaxSpringLength") == 0)
            mSpringCutoffLength = atof(fields[param_value_index].c_str());
        
        cout << "Parameter: " << fields[param_name_index] << ", has value: " << fields[param_value_index] << endl;
	}
}


void BasicProcess::loadMolecule() {
    
    CSVRead molecule_file;
    
    /*=================================================
            Identify number of atoms in the molecule
     */
    if(!molecule_file.OpenCSVFile(mMoleculeFile))
        throw "Unable to open molecular file";
    
    vector<string> molecule_columns;
	int x_index   = static_cast<int>(molecule_columns.size()); molecule_columns.push_back("X" );
    
    if(!molecule_file.ReadCSVHeader(molecule_columns)) {
		throw "Molecular file didn't have correct column headers";
	}
    
    int max_index = 0;
	while(true) {
		vector<string> fields = molecule_file.ReadCSVRecord();
		if(fields.size() == 0) {
			break;
		}
        
        ++max_index;
	}
    
    if(max_index <= 0)
        throw "No atoms";
    
    cout << mMoleculeFile << " contains a total of " << max_index << " atoms." << endl;
    mMolecule.AllocateNodes(max_index);
    /*==================================================
     */
    
    
    mMolecule.PrintLimits();
    /*=================================================
            Save the atoms as nodes to the molecule
     */
    if(!molecule_file.OpenCSVFile(mMoleculeFile))
        throw "Unable to open molecular file";
    
    int y_index       = static_cast<int>(molecule_columns.size()); molecule_columns.push_back("Y");
    int z_index       = static_cast<int>(molecule_columns.size()); molecule_columns.push_back("Z");
    
	if(!molecule_file.ReadCSVHeader(molecule_columns)) {
		throw "Molecular file didn't have correct column headers";
	}
     
    while(true) {
		vector<string> fields = molecule_file.ReadCSVRecord();
		if(fields.size() == 0) {
			break;
		}
        
        mMolecule.AddNode(
            atof(fields[x_index].c_str()),
            atof(fields[y_index].c_str()),
            atof(fields[z_index].c_str()) );
	}
    
    mMolecule.PrintLimits();
}