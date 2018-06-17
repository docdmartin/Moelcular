#include <map>
#include <limits>
#include <iostream>
#include <cmath>

#include "Molecule.h"

using namespace std;

struct coord {
    int ix, iy, iz;
    
    bool operator==(const coord &o) const {
        return (ix == o.ix && iy == o.iy && iz == o.iz);
    }
    
    bool operator<(const coord &o) const {
        return (ix < o.ix || (ix == o.ix && iy < o.iy) || (ix == o.ix && iy == o.iy && iz < o.iz));
    }
};



Molecule::Molecule() :
  mMinX( numeric_limits<double>::max() ), mMaxX( -numeric_limits<double>::max() ),
  mMinY( numeric_limits<double>::max() ), mMaxY( -numeric_limits<double>::max() ),
  mMinZ( numeric_limits<double>::max() ), mMaxZ( -numeric_limits<double>::max() )
{
    
}

Molecule::~Molecule(){
    
}

void Molecule::AddNode(double x, double y, double z) {
    if(mMinX > x)
        mMinX = x;
    if(mMaxX < x)
        mMaxX = x;
    
    if(mMinY > y)
        mMinY = y;
    if(mMaxY < y)
        mMaxY = y;
    
    if(mMinZ > z)
        mMinZ = z;
    if(mMaxZ < z)
        mMaxZ = z;
    
    mNodes.push_back(Node(static_cast<int>(mNodes.size()), x, y, z));
}

void Molecule::IdentifyContacts( double cutoff ) {
    
    
    map<coord, vector<int> > voxels;
    map<coord, vector<int> >::iterator vit;
    
    vector<Node>::iterator it;
    for(it = mNodes.begin(); it != mNodes.end(); ++it) {
        int ix = static_cast<int>( floor( (it->GetX() - mMinX) / cutoff ) );
        int iy = static_cast<int>( floor( (it->GetY() - mMinY) / cutoff ) );
        int iz = static_cast<int>( floor( (it->GetZ() - mMinZ) / cutoff ) );
        
        vector<int> tmp_node_ids;
        
        coord key = {ix, iy, iz};
        vit = voxels.find( key );
        
        if(vit != voxels.end()) {
            tmp_node_ids = vit->second;
        }
        tmp_node_ids.push_back( it->GetID() );
        
        voxels.insert( pair<coord, vector<int> >(key, tmp_node_ids) );
    }
    
    vector<int>::iterator nit;
    vector<int>::iterator niter;
    for(vit = voxels.begin(); vit != voxels.end(); ++vit) {
        for(nit = vit->second.begin(); nit != vit->second.end(); ++nit) {
            // First check nodes in same voxel
            for(niter = nit + 1; niter != vit->second.end(); ++niter) {
                // check radius and construct connection if needed
            }
            
            // check 7 neighboring voxels (only increasing voxels)
            //   x+1, y, z
            //   x, y+1, z
            //   x+1, y+1, z
            //   x, y, z+1
            //   x+1, y, z+1
            //   x, y+1, z+1
            //   x+1, y+1, z+1
            //   build new key
            //   get node vector from voxels
            //   test radius for connections
        }
    }
}


void Molecule::PrintLimits() {
cout << "X: (" << mMinX << ", " << mMaxX << ")" << endl;
cout << "Y: (" << mMinY << ", " << mMaxY << ")" << endl;
cout << "Z: (" << mMinZ << ", " << mMaxZ << ")" << endl;
}