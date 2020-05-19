#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

class Configuration {

public:

  Configuration(){}
  ~Configuration(){}

  bool   mNoncovalentRangeDep = false;

  int    mNumSegments         =   50;//100;
  double mMinFreq             =   0.1;//  0.001;
  double mMaxFreq             = 100.;//10000.0;
  double mChargeConst         =   1.;//332.063182082;
  double mCovalentStrength    =   10.;//100.;
  double mNoncovalentStrength =     0.9;//1.;
  double mNoncovalentRange    =    2.;//15.;

  string mNodeCentroid  = "CA";

private:

};

#endif
