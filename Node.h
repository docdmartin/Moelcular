#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include "Common.h"
#include "Atom.h"

using namespace std;

class Node {

public:

  Node(){}
  ~Node(){}

  void SetPosition( double x, double y, double z ){
    mPosition[0] = x;
    mPosition[1] = y;
    mPosition[2] = z;
  }
  void SetPositionCarbonAlpha(){
    for( Atom* it = mBegin; it != mEnd; ++it ){
      if( (*it).GetName() == ElementType::C && (*it).GetGreek() == 'A' ){
        double* pos = (*it).GetPosition();
        mPosition[0] = *(pos);
        mPosition[1] = *(pos+1);
        mPosition[2] = *(pos+2);
        return;
      }
    }

    Print();
    throw string("Alpha Carbon not found");
  }
  double* GetPosition() { return mPosition; }

  void SetName( string& name ) { mName = name; }
  string& GetName() { return mName; }

  void SetId( string& name ) { mId = name; }
  string& GetId() { return mId; }

  // Assumes atoms are processed in order
  void AddAtom( Atom* it ){
    if( mBegin == nullptr )
      mBegin = it;

    mEnd = it;
  }
  Atom* GetBegin() { return mBegin; }
  Atom* GetEnd  () { return mEnd  ; }

  void SetElectricPotential( double* E, double* pos )
  {
    E[0] = 0.; E[1] = 0.; E[2] = 0.;
    double dp[3], c;
    double* p;

    dp[0] = pos[0] - mPosition[0];
    dp[1] = pos[1] - mPosition[1];
    dp[2] = pos[2] - mPosition[2];
    c = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
    if( c < 1e-12 )
      return;

    for( Atom* it = mBegin; it <= mEnd; ++it ){
      p = (*it).GetPosition();
      dp[0] = pos[0] - p[0];
      dp[1] = pos[1] - p[1];
      dp[2] = pos[2] - p[2];

      c = std::sqrt( dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2] );
      if( c < 1e-15 ) // Probe was placed at location of an atom, atom is being ignored in electric potential
       continue;

      c = (*it).GetCharge() / ( c * c * c )   * 332.063182082;
      E[0] += dp[0] * c;
      E[1] += dp[1] * c;
      E[2] += dp[2] * c;
    }

  //  cout << "Force: " << E[0] << ", " << E[1] << ", " << E[2] << endl;
  }


  void Print(){
    cout << mName << " at location: (" << mPosition[0] << ", " << mPosition[1] << ", " << mPosition[2] << ")" << endl;
    if( mBegin == nullptr || mEnd == nullptr )
      return;

    for( Atom* it = mBegin; it != mEnd; ++it ){
      cout << "\t";
      (*it).Print();
    }
  }

private:

  double mPosition[3];
  Atom*  mBegin = nullptr;
  Atom*  mEnd   = nullptr;
  string mName;
  string mId;
};

#endif
