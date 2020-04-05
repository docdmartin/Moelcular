#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include "Common.h"

using namespace std;

class Atom {

public:

  Atom(){}
  ~Atom(){}

  void SetPosition( double x, double y, double z );
  double* GetPosition() { return mPosition; }

  void SetCharge( double q ) { mCharge = q; }
  double GetCharge() { return mCharge; }

  void SetName( string& name );
  ElementType GetName() { return mName; }

  char GetGreek() { return mGreek; }

  void Print(){
    cout << EnumToString( mName ) << " with charge = " << mCharge
         << " at location: (" << mPosition[0] << ", " << mPosition[1] << ", " << mPosition[2] << ")" << endl;
  }

private:

  double      mPosition[3];
  double      mCharge;
  ElementType mName;
  char        mGreek;

};

#endif
