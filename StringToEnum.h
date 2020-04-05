#include<string>
#include<map>
#include<algorithm>

void StringToEnum( std::string& str, ElementType& value ){
  std::map< std::string, ElementType >::iterator it = gElementTypeMap.begin();
  for(; it != gElementTypeMap.end(); ++it){
    std::string tmp = it->first;
    if( tmp == str ){
      value = it->second;
      return;
    }
  }

  value = ElementSize;
}
