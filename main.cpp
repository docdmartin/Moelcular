#include <string>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <exception>

#include <map>
#include <vector>

#include "ElasticNetworkModel.h"
#include "util/CSVRead.h"

using namespace std;

void loadConfigurationFile(string configurationFile);
void loadNetworkModel     (string proteinFolderName);
int  countAminoAcid       (string input_file, string chain_name);
void loadAminoAcid        (string input_file, string chain_name);

ElasticNetworkModel networkModel;
map<string, string> configurationMap;

int main(int argc,char **argv)
{
    if(argc == 1){
        cout << "Configuration file must be included on command line." << endl;
        cout << "> " << argv[0] << " pathname/filename" << endl;

        return 1;
    }

    string filename = argv[1];

    string executable_file = argv[0];
    size_t botDirPos = executable_file.find_last_of("/");
    string current_working_dir = executable_file.substr(0, botDirPos+1);
    cout << "Local directory: " << current_working_dir << endl;


    loadConfigurationFile( filename );


    map<string, string>::iterator it = configurationMap.find("MaxSpringLength");
    if(it == configurationMap.end()) {
        cout << "Configuration file did not include: MaxSpringLength" << endl;
        return 0;
    }
    networkModel.SetMaxSpringLength( atof( it->second.c_str() ) );


    it = configurationMap.find("NodeType");
    if(it == configurationMap.end()) {
        cout << "Configuration file did not include: NodeType" << endl;
        return 0;
    }
    networkModel.SetNodeType( it->second );


    it = configurationMap.find("ProteinFolderName");
    if(it == configurationMap.end()) {
        cout << "Configuration file did not include: ProteinFolderName" << endl;
        return 0;
    }
    loadNetworkModel( it->second );

    networkModel.ConfigureModel();


    try {
        networkModel.IdentifyContacts();
        networkModel.ConstructLinearResponse();

        int ref_node_index = networkModel.AddReferencePoint( 0.0, 0.0, 0.0 );

        vector< pair<double, double> > complex_potential1 = networkModel.SingleModeFrequencyResponse(ref_node_index, 1.0);
        cout << "Electric potential result = " << complex_potential1[0].first << " + " << complex_potential1[0].second << " i " << endl;
        cout << "Dipole potential result = " << complex_potential1[1].first << " + " << complex_potential1[1].second << " i " << endl;

        vector< pair<double, double> > complex_potential2 = networkModel.DualModeFrequencyResponse(ref_node_index, 1.0, .3, 0.5);
        cout << "Electric potential result = " << complex_potential2[0].first << " + " << complex_potential2[0].second << " i " << endl;
        cout << "Dipole potential result = " << complex_potential2[1].first << " + " << complex_potential2[1].second << " i " << endl;
    }
    catch (const char* msg) {
        cerr << msg << endl;
    }

	  return 1;
}


void loadConfigurationFile(string configurationFile){

    CSVRead config_file;

    if(!config_file.OpenCSVFile(configurationFile))
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

        configurationMap[fields[param_name_index]] = fields[param_value_index];

        cout << "Parameter: " << fields[param_name_index] << ", has value: " << fields[param_value_index] << endl;
	  }
}



void loadNetworkModel(string proteinFolderName) {

    DIR           *dir;
    struct dirent *ent;
    if ((dir = opendir (proteinFolderName.c_str())) != NULL) {

        size_t botDirPos   = proteinFolderName.find_last_of("/");
        string proteinName = proteinFolderName.substr(botDirPos+1, proteinFolderName.size() - botDirPos);

        cout << "Protein name: " << proteinName << endl;

        /* print all the files and directories within directory */
        int amino_cnt = 0;
        while ((ent = readdir (dir)) != NULL) {
            string tmpFileName = ent->d_name;

            if (tmpFileName.compare(0, proteinName.size(), proteinName))
                continue;

            string chain = tmpFileName.substr(proteinName.size()+1, tmpFileName.size() - proteinName.size() - 5);

            amino_cnt += countAminoAcid(proteinFolderName + "/" + tmpFileName, chain);
        }
        closedir (dir);


        if(amino_cnt <= 0)
            throw "No nodes";

        cout << proteinFolderName << " contains a total of " << amino_cnt << " alpha carbons." << endl;
        networkModel.AllocateNodes( amino_cnt );


        /* print all the files and directories within directory */
        dir = opendir (proteinFolderName.c_str());
        while ((ent = readdir (dir)) != NULL) {
            string tmpFileName = ent->d_name;

            if (tmpFileName.compare(0, proteinName.size(), proteinName))
                continue;

            string chain = tmpFileName.substr(proteinName.size()+1, tmpFileName.size() - proteinName.size() - 5);

            loadAminoAcid(proteinFolderName + "/" + tmpFileName, chain);
        }
        closedir (dir);

    }
}

int countAminoAcid(string input_file, string chain_name) {

    CSVRead network_file;

    if(!network_file.OpenCSVFile(input_file))
        throw "Unable to open network file";

    vector<string> network_columns;
    int name_index   = static_cast<int>(network_columns.size()); network_columns.push_back("ATOMNAME" );

    if(!network_file.ReadCSVHeader(network_columns)) {
        throw "Network file didn't have correct column headers";
    }

    int cnt = 0;
    while(true) {
        vector<string> fields = network_file.ReadCSVRecord();
        if(fields.size() == 0) {
            break;
        }

        // If atom is not alpha carbon then go to next line
        if(fields[name_index].compare("CA") != 0)
            continue;

        ++cnt;
    }

    return cnt;
}

void loadAminoAcid(string input_file, string chain_name) {

    cout << "Loading " << chain_name << " from file: " << input_file << endl;

    CSVRead network_file;

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
            node_index          = networkModel.AddNode(resid, resname, atomid, atomname, x, y, z, q);

            if(node_index < 0)
                throw "Process terminated due to error";

            if(prev_node_index >= 0)
                networkModel.CreateBackboneConnection(prev_node_index, node_index);
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
        node_index          = networkModel.AddNode(resid, resname, atomid, atomname, x, y, z, q);

        if(prev_node_index >= 0)
            networkModel.CreateBackboneConnection(prev_node_index, node_index);
    }
}
