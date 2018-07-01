#ifndef ____BasicProcess__
#define ____BasicProcess__

#include <iostream>
#include <string>


#include "network_model/Network.h"

using namespace std;

class BasicProcess{
public:
    BasicProcess(string);
    ~BasicProcess();

private:
    void loadConfigurationFile();
    void loadNetworkModel();

    string  mConfigurationFile;
    string  mInputFile;

    Network mNetworkModel;

};

#endif
