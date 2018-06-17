#include "Node.h"

Node::Node(int id, double x, double y, double z)
  : mNodeID(id)
{
    mVoxelIndex[0] = -1; mVoxelIndex[1] = -1; mVoxelIndex[2] = -1;
    mPos       [0] =  x; mPos       [1] =  y; mPos       [2] =  z;
}

Node::~Node(){
    
}

void Node::SetPosition(double x, double y, double z) {
    mPos[0] = x;
    mPos[1] = y;
    mPos[2] = z;
}

void Node::SetIndex(int x, int y, int z) {
    mVoxelIndex[0] = x;
    mVoxelIndex[1] = y;
    mVoxelIndex[2] = z;
}