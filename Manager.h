#ifndef MANAGER_H
#define MANAGER_H

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "Configuration.h"
#include "Node.h"
#include "Atom.h"

using namespace std;

class Manager {

public:

  Manager(){}
  ~Manager(){}

  void CreateConfiguration( int id ) {
    if( mConfigMap.empty() )
      mId = id;
    mConfigMap.insert( pair<int, Configuration>(id, Configuration()) );
  }

  void SetNodeCount( int cnt ) { mNodeCount = cnt; mNodeArray.reserve(mNodeCount); }
  void SetAtomCount( int cnt ) { mAtomCount = cnt; mAtomArray.reserve(mAtomCount); }

  void SetParameter( int id, char* param_name, char* param_value ){
    string param_str = param_name;

    if( param_str == "Min Freq" ){
      mConfigMap[id].mMinFreq = static_cast<double>( atof(param_value) );
    }
    else if( param_str == "Max Freq" ){
      mConfigMap[id].mMaxFreq = static_cast<double>( atof(param_value) );
    }
    else if( param_str == "Num Freq Segments" ){
      mConfigMap[id].mNumSegments = static_cast<int>( atoi(param_value) );
    }
    else if( param_str == "Charge Const" ){
      mConfigMap[id].mChargeConst = static_cast<double>( atof(param_value) );
    }
    else if( param_str == "Covalent Strength" ){
      mConfigMap[id].mCovalentStrength = static_cast<double>( atof(param_value) );
    }
    else if( param_str == "Non-covalent Strength" ){
      mConfigMap[id].mNoncovalentStrength = static_cast<double>( atof(param_value) );
    }
    else if( param_str == "Non-covalent Max Range" ){
      mConfigMap[id].mNoncovalentRange = static_cast<double>( atof(param_value) );
    }
    else if( param_str == "Non-covalent Range Dep" ){
      if( static_cast<int>( atoi(param_value) ) == 0 )
        mConfigMap[id].mNoncovalentRangeDep = false;
      else
        mConfigMap[id].mNoncovalentRangeDep = true;
    }
    else if( param_str == "Node Centroid Type" ){
      mConfigMap[id].mNodeCentroid = param_value;
    }
  }

  void AddAtom( string& segid, string& resid, string& atomname, double* pos, double mass, double charge ){
    mAtomArray.push_back( Atom() );
    mAtomArray.back().SetPosition(pos[0], pos[1], pos[2]);
    mAtomArray.back().SetCharge( charge * GetChargeConst() );
    mAtomArray.back().SetName( atomname );
    mAtomArray.back().SetMass(mass);
    mAtomArray.back().SetSegment(segid);
    mAtomArray.back().SetResId( resid );
  }

  vector<double> GetFrequency(){
    vector<double> freq;
    freq.reserve( mConfigMap[mId].mNumSegments+1 );

    double log_min = log( mConfigMap[mId].mMinFreq );
    double log_max = log( mConfigMap[mId].mMaxFreq );
    double mult    = exp( (log_max - log_min)/static_cast<double>(mConfigMap[mId].mNumSegments) );

    freq.push_back( mConfigMap[mId].mMaxFreq );
    double tmp = mConfigMap[mId].mMaxFreq;
    for( int cnt = 0; cnt < mConfigMap[mId].mNumSegments; ++cnt ){
      tmp /= mult;
      freq.push_back( tmp );
    }

    return freq;
  }

  double GetChargeConst() { return mConfigMap[mId].mChargeConst; }
  double GetCovalentStrength() { return mConfigMap[mId].mCovalentStrength; }
  double GetNoncovalentStrength() { return mConfigMap[mId].mNoncovalentStrength; }
  double GetNoncovalentRange() { return mConfigMap[mId].mNoncovalentRange; }
  bool   IsBondStrengthRangeDep() { return mConfigMap[mId].mNoncovalentRangeDep; }
  string CentroidType() { return mConfigMap[mId].mNodeCentroid; }

  int GetNodeCount() { return mNodeCount; }
  int GetAtomCount() { return mAtomCount; }

  void BuildNodes(){
    mNodeArray.push_back( Node() );
    string prev_segid = mAtomArray[0].GetSegment();
    string prev_resid = mAtomArray[0].GetResId  ();
    for( int cnt = 0; cnt < mAtomCount; ++cnt ){
      if( prev_segid != mAtomArray[cnt].GetSegment() || prev_resid != mAtomArray[cnt].GetResId() ){
        mNodeArray.push_back( Node() );

        if( prev_segid == mAtomArray[cnt].GetSegment() ){
          pair<size_t, size_t> link(mNodeArray.size()-2, mNodeArray.size()-1);
          mLinkMap.insert( pair<pair<size_t, size_t>, double>(link, GetCovalentStrength()) );
        }

        prev_segid = mAtomArray[cnt].GetSegment();
        prev_resid = mAtomArray[cnt].GetResId  ();
      }
      mNodeArray.back().AddAtom( &(mAtomArray[cnt]) );
    }

    for( auto& node : mNodeArray )
      node.SetPositionCarbonAlpha();


    double dx, dy, dz;
    double r = GetNoncovalentRange();
    double rr = r*r;
    for( int node1 = 0; node1 < mNodeCount-1; ++node1 ){
      double* node1_pos = mNodeArray[node1].GetPosition();

      for( int node2 = node1+1; node2 < mNodeCount; ++node2 ){
        if( node2 == node1+1 && mNodeArray[node2].GetSegment() == mNodeArray[node1].GetSegment() )
          continue;

        double* node2_pos = mNodeArray[node2].GetPosition();

        dx = *(node1_pos  ) - *(node2_pos  );
        dy = *(node1_pos+1) - *(node2_pos+1);
        dz = *(node1_pos+2) - *(node2_pos+2);
        if( fabs(dx) > r || fabs(dy) > r || fabs(dz) > r )
          continue;

        if( dx*dx + dy*dy + dz*dz < rr ){
          pair<size_t, size_t> link(node1, node2);
          mLinkMap.insert( pair<pair<size_t, size_t>, double>(link, GetNoncovalentStrength()) );
        }
      }
    }
  }

  void AddProbe( int id, double* pos ){
    vector<double> loc;
    loc.push_back( pos[0] );
    loc.push_back( pos[1] );
    loc.push_back( pos[2] );
    mProbe.insert( pair<int, vector<double>>(id, loc) );
  }

  void AddKernel( int id, double weight, double freq_mult ){
    (void) id;
    mKernel.push_back( pair<double, double>(weight, freq_mult) );
  }

  vector<Node>& GetNodeArray() { return mNodeArray; }
  map<pair<size_t, size_t>, double>& GetLinkMap() { return mLinkMap; }
  map< int, vector<double> >& GetProbe() { return mProbe; }
  vector< pair<double, double> >& GetKernel() {
    if( mKernel.size() == 0 ){
      mKernel.push_back( pair<double, double>(1., 1.) );
    }
    return mKernel;
  }

private:

  int mId;
  int mNodeCount = 0;
  int mAtomCount = 0;

  vector<Node> mNodeArray;
  vector<Atom> mAtomArray;

  map< int                 , Configuration        > mConfigMap;
  map< pair<size_t, size_t>, double               > mLinkMap;
  map< int                 , vector<double>       > mProbe;
  vector< pair<double, double>                    > mKernel;

};

#endif
