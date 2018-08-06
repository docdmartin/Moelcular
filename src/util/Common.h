#ifndef ____Common__
#define ____Common__

#include <string>
#include <set>
#include <map>
#include <vector>
#include "util/CommonType.h"

using namespace std;

class Common{
public:
    Common();
    ~Common();

    void SetNodeType(CommonType::NodeType n) { mNodeType = n; }
    CommonType::NodeType GetNodeType() { return mNodeType; }

    void SetVariables();

    struct ElementData{
      string abbreviation;
      string name;
      double mass;
    };


    map<CommonType::ElementType, ElementData> mPeriodicTable;
    map<CommonType::NodeType, string>         mNodeTypeDef;

private:
    CommonType::NodeType                      mNodeType;
};

#endif
