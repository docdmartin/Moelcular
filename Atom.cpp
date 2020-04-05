#include "Atom.h"

using namespace std;

void Atom::SetPosition( double x, double y, double z ){
  mPosition[0] = x;
  mPosition[1] = y;
  mPosition[2] = z;
}

void Atom::SetName( string& name ){

  char first_letter = name.at(0);
  if( name.size() > 1 ){
    mGreek = name.at(1);
  }
  else{
    mGreek = '.';
  }

  string tmp(1, first_letter);
  StringToEnum( tmp, mName );
}
