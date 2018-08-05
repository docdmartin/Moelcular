#ifndef ____Common__
#define ____Common__


namespace CommonEnum {
  /*
    ConnectionType must has NO_CONNECTION specified as the first in a list of enum
    and ALL_CONNECTION as the last. This is important since we iterate over the enum list at times
  */
  enum ConnectionType
  {
    NO_CONNECTION,
    SPRING_LEVEL_1,
    SPRING_LEVEL_2,
    ALL_CONNECTION
  };

  enum NodeType{
    ALPHA_CARBON,
    AUXILIARY
  };
}

#endif
