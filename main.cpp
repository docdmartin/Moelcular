#include <string>
#include <iostream>
#include <exception>

#include "src/process/BasicProcess.h"

using namespace std;

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

    try {
        BasicProcess process( filename, current_working_dir );
    }
    catch (const char* msg) {
        cerr << msg << endl;
    }

	return 1;
}
