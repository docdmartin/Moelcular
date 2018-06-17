#ifndef ____BasicProcess__
#define ____BasicProcess__

#include <iostream>
#include <string>


#include "../../src/model/Molecule.h"

using namespace std;

class BasicProcess{
public:
    BasicProcess(string);
    ~BasicProcess();
    
private:
    void loadConfigurationFile();
    void loadMolecule();
    
    string mConfigurationFile;
    string mMoleculeFile;
    
    Molecule mMolecule;
    
    double mSpringCutoffLength;
    
};

#endif
