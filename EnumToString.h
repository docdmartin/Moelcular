#include<string>
#include<map>

std::string EnumToString( ElementType value ){
  std::map< std::string, ElementType >::iterator it = gElementTypeMap.begin();
  for(; it != gElementTypeMap.end(); ++it){
    if( it->second == value )
      return it->first;
  }

  return std::string("Unknown");
}
