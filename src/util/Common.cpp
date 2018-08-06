#include "util/Common.h"

Common::Common(){
  SetVariables();
}

Common::~Common(){
}

void Common::SetVariables() {
      mNodeType = CommonType::NodeType::UNDEFINED;

      mNodeTypeDef.insert( pair<CommonType::NodeType, string>(CommonType::NodeType::UNDEFINED         , "Undefined") );
      mNodeTypeDef.insert( pair<CommonType::NodeType, string>(CommonType::NodeType::ALPHA_CARBON      , "Alpha Carbon") );
      mNodeTypeDef.insert( pair<CommonType::NodeType, string>(CommonType::NodeType::MASS_WEIGHTED_MEAN, "Mass Weighted Mean") );

      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::H , ElementData{"H" , "Hydrogen"     ,   1.008      }) );
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::He, ElementData{"He", "Helium"       ,   4.002602   }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Li, ElementData{"Li", "Lithium"      ,   6.94       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Be, ElementData{"Be", "Beryllium"    ,   9.0121831  }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::B , ElementData{"B" , "Boron"        ,  10.81       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::C , ElementData{"C" , "Carbon"       ,  12.011      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::N , ElementData{"N" , "Nitrogen"     ,  14.007      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::O , ElementData{"O" , "Oxygen"       ,  15.999      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::F , ElementData{"F" , "Fluorine"     ,  18.998403163}));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ne, ElementData{"Ne", "Neon"         ,  20.1797     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Na, ElementData{"Na", "Sodium"       ,  22.98976928 }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Mg, ElementData{"Mg", "Magnesium"    ,  24.305      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Al, ElementData{"Al", "Aluminium"    ,  26.9815385  }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Si, ElementData{"Si", "Silicon"      ,  28.085      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::P , ElementData{"P" , "Phosphorus"   ,  30.973761998}));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::S , ElementData{"S" , "Sulfur"       ,  32.06       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cl, ElementData{"Cl", "Chlorine"     ,  35.45       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ar, ElementData{"Ar", "Argon"        ,  39.948      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::K , ElementData{"K" , "Potassium"    ,  39.0983     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ca, ElementData{"Ca", "Calcium"      ,  40.078      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Sc, ElementData{"Sc", "Scandium"     ,  44.955908   }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ti, ElementData{"Ti", "Titanium"     ,  47.867      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::V , ElementData{"V" , "Vanadium"     ,  50.9415     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cr, ElementData{"Cr", "Chromium"     ,  51.9961     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Mn, ElementData{"Mn", "Manganese"    ,  54.938044   }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Fe, ElementData{"Fe", "Iron"         ,  55.845      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Co, ElementData{"Co", "Cobalt"       ,  58.933194   }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ni, ElementData{"Ni", "Nickel"       ,  58.6934     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cu, ElementData{"Cu", "Copper"       ,  63.546      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Zn, ElementData{"Zn", "Zinc"         ,  65.38       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ga, ElementData{"Ga", "Gallium"      ,  69.723      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ge, ElementData{"Ge", "Germanium"    ,  72.630      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::As, ElementData{"As", "Arsenic"      ,  74.921595   }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Se, ElementData{"Se", "Selenium"     ,  78.971      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Br, ElementData{"Br", "Bromine"      ,  79.904      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Kr, ElementData{"Kr", "Krypton"      ,  83.798      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Rb, ElementData{"Rb", "Rubidium"     ,  85.4678     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Sr, ElementData{"Sr", "Strontium"    ,  87.62       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Y , ElementData{"Y" , "Yttrium"      ,  88.90584    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Zr, ElementData{"Zr", "Zirconium"    ,  91.224      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Nb, ElementData{"Nb", "Niobium"      ,  92.90637    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Mo, ElementData{"Mo", "Molybdenum"   ,  95.95       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Tc, ElementData{"Tc", "Technetium"   ,  98          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ru, ElementData{"Ru", "Ruthenium"    , 101.07       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Rh, ElementData{"Rh", "Rhodium"      , 102.90550    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pd, ElementData{"Pd", "Palladium"    , 106.42       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ag, ElementData{"Ag", "Silver"       , 107.8682     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cd, ElementData{"Cd", "Cadmium"      , 112.414      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::In, ElementData{"In", "Indium"       , 114.818      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Sn, ElementData{"Sn", "Tin"          , 118.710      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Sb, ElementData{"Sb", "Antimony"     , 121.760      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Te, ElementData{"Te", "Tellurium"    , 127.60       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::I , ElementData{"I" , "Iodine"       , 126.90447    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Xe, ElementData{"Xe", "Xenon"        , 131.293      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cs, ElementData{"Dc", "Caesium"      , 132.90545196 }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ba, ElementData{"Ba", "Barium"       , 137.327      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::La, ElementData{"La", "Lanthanum"    , 138.90547    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ce, ElementData{"Ce", "Cerium"       , 140.116      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pr, ElementData{"Pr", "Praseodymium" , 140.90766    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Nd, ElementData{"Nd", "Neodymium"    , 144.242      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pm, ElementData{"Pm", "Promethium"   , 145          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Sm, ElementData{"Sm", "Samarium"     , 150.36       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Eu, ElementData{"Eu", "Europium"     , 151.964      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Gd, ElementData{"Gd", "Gadolinium"   , 157.25       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Tb, ElementData{"Tb", "Terbium"      , 158.92535    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Dy, ElementData{"Dy", "Dysprosium"   , 162.500      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ho, ElementData{"Ho", "Holmium"      , 164.93033    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Er, ElementData{"Er", "Erbium"       , 167.259      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Tm, ElementData{"Tm", "Thulium"      , 168.93422    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Yb, ElementData{"Yb", "Ytterbium"    , 173.045      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Lu, ElementData{"Lu", "Lutetium"     , 174.9668     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Hf, ElementData{"Hf", "Hafnium"      , 178.49       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ta, ElementData{"Ta", "Tantalum"     , 180.94788    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::W , ElementData{"W" , "Tungsten"     , 183.84       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Re, ElementData{"Re", "Rhenium"      , 186.207      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Os, ElementData{"Os", "Osmium"       , 190.23       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ir, ElementData{"Ir", "Iridium"      , 192.217      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pt, ElementData{"Pt", "Platinum"     , 195.084      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Au, ElementData{"Au", "Gold"         , 196.966569   }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Hg, ElementData{"Hg", "Mercury"      , 200.592      }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Tl, ElementData{"Tl", "Thallium"     , 204.38       }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pb, ElementData{"Pb", "Lead"         , 207.2        }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Bi, ElementData{"Bi", "Bismuth"      , 208.98040    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Po, ElementData{"Po", "Polonium"     , 209          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::At, ElementData{"At", "Astatine"     , 210          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Rn, ElementData{"Rn", "Radon"        , 222          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Fr, ElementData{"Fr", "Francium"     , 223          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ra, ElementData{"Ra", "Radium"       , 226          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ac, ElementData{"Ac", "Actinium"     , 227          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Th, ElementData{"Th", "Thorium"      , 232.0377     }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pa, ElementData{"Pa", "Protactinium" , 231.03588    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::U , ElementData{"U" , "Uranium"      , 238.02891    }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Np, ElementData{"Np", "Neptunium"    , 237          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Pu, ElementData{"Pu", "Plutonium"    , 244          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Am, ElementData{"Am", "Americium"    , 243          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cm, ElementData{"Cm", "Curium"       , 247          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Bk, ElementData{"Bk", "Berkelium"    , 247          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cf, ElementData{"Cf", "Californium"  , 251          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Es, ElementData{"Es", "Einsteinium"  , 252          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Fm, ElementData{"Fm", "Fermium"      , 257          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Md, ElementData{"Md", "Mendelevium"  , 258          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::No, ElementData{"No", "Nobelium"     , 259          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Lr, ElementData{"Lr", "Lawrencium"   , 266          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Rf, ElementData{"Rf", "Rutherfordium", 267          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Db, ElementData{"Db", "Dubnium"      , 268          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Sg, ElementData{"Sg", "Seaborgium"   , 269          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Bh, ElementData{"Bh", "Bohrium"      , 270          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Hs, ElementData{"Hs", "Hassium"      , 277          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Mt, ElementData{"Mt", "Meitnerium"   , 278          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ds, ElementData{"Ds", "Darmstadtium" , 281          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Rg, ElementData{"Rg", "Roentgenium"  , 282          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Cn, ElementData{"Cn", "Copernicium"  , 285          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Nh, ElementData{"Nh", "Nihonium"     , 286          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Fl, ElementData{"Fl", "Flerovium"    , 289          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Mc, ElementData{"Mc", "Moscovium"    , 290          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Lv, ElementData{"Lv", "Livermorium"  , 293          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Ts, ElementData{"Ts", "Tennessine"   , 294          }));
      mPeriodicTable.insert( pair<CommonType::ElementType, ElementData>(CommonType::ElementType::Og, ElementData{"Og", "Oganesson"    , 294          }));

}
