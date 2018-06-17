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
    
    try {
        BasicProcess process( filename );
    }
    catch (const char* msg) {
        cerr << msg << endl;
    }
    
	return 1;
}

