#include <string>
#include <iostream>
#include <exception>

#include <Eigen/Dense>

#include "src/process/BasicProcess.h"

using namespace std;

int main(int argc,char **argv)
{
    Eigen::MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    cout << m << endl;


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
