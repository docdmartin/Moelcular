#ifndef ____BasicProcess__
#define ____BasicProcess__

#include <iostream>
#include <string>


#include "network_model/Network.h"

using namespace std;

class BasicProcess{
public:
    BasicProcess(string, string);
    ~BasicProcess();

private:
    void loadConfigurationFile();
    void loadNetworkModel();
    void loadAminoAcid(string, string);
    void countAminoAcid(string, string);

    int     mMaxIndex;
    string  mConfigurationFile;
    string  mPathName;

    string  mInputFolder;

    Network mNetworkModel;

    Common  mParameter;
};

#endif
