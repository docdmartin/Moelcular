#include "util/Common.h"

using namespace CommonType;

Common::Common(){
  SetVariables();
}

Common::~Common(){
}

void Common::SetVariables() {
      mNodeType = NodeType::UNDEFINED;

      mConnectionStrength[CommonType::NO_CONNECTION ] =  0.0;
      mConnectionStrength[CommonType::SPRING_LEVEL_1] = 10.0;
      mConnectionStrength[CommonType::SPRING_LEVEL_2] =  1.0;
      mConnectionStrength[CommonType::ALL_CONNECTION] =  1.0;

      mNodeTypeDef[NodeType::UNDEFINED         ] = "Undefined";
      mNodeTypeDef[NodeType::ALPHA_CARBON      ] = "Alpha Carbon";
      mNodeTypeDef[NodeType::MASS_WEIGHTED_MEAN] = "Mass Weighted Mean";

      mPeriodicTable[H ] = ElementData{"H" , "Hydrogen"     ,   1.008      };
      mPeriodicTable[He] = ElementData{"He", "Helium"       ,   4.002602   };
      mPeriodicTable[Li] = ElementData{"Li", "Lithium"      ,   6.94       };
      mPeriodicTable[Be] = ElementData{"Be", "Beryllium"    ,   9.0121831  };
      mPeriodicTable[B ] = ElementData{"B" , "Boron"        ,  10.81       };
      mPeriodicTable[C ] = ElementData{"C" , "Carbon"       ,  12.011      };
      mPeriodicTable[N ] = ElementData{"N" , "Nitrogen"     ,  14.007      };
      mPeriodicTable[O ] = ElementData{"O" , "Oxygen"       ,  15.999      };
      mPeriodicTable[F ] = ElementData{"F" , "Fluorine"     ,  18.998403163};
      mPeriodicTable[Ne] = ElementData{"Ne", "Neon"         ,  20.1797     };
      mPeriodicTable[Na] = ElementData{"Na", "Sodium"       ,  22.98976928 };
      mPeriodicTable[Mg] = ElementData{"Mg", "Magnesium"    ,  24.305      };
      mPeriodicTable[Al] = ElementData{"Al", "Aluminium"    ,  26.9815385  };
      mPeriodicTable[Si] = ElementData{"Si", "Silicon"      ,  28.085      };
      mPeriodicTable[P ] = ElementData{"P" , "Phosphorus"   ,  30.973761998};
      mPeriodicTable[S ] = ElementData{"S" , "Sulfur"       ,  32.06       };
      mPeriodicTable[Cl] = ElementData{"Cl", "Chlorine"     ,  35.45       };
      mPeriodicTable[Ar] = ElementData{"Ar", "Argon"        ,  39.948      };
      mPeriodicTable[K ] = ElementData{"K" , "Potassium"    ,  39.0983     };
      mPeriodicTable[Ca] = ElementData{"Ca", "Calcium"      ,  40.078      };
      mPeriodicTable[Sc] = ElementData{"Sc", "Scandium"     ,  44.955908   };
      mPeriodicTable[Ti] = ElementData{"Ti", "Titanium"     ,  47.867      };
      mPeriodicTable[V ] = ElementData{"V" , "Vanadium"     ,  50.9415     };
      mPeriodicTable[Cr] = ElementData{"Cr", "Chromium"     ,  51.9961     };
      mPeriodicTable[Mn] = ElementData{"Mn", "Manganese"    ,  54.938044   };
      mPeriodicTable[Fe] = ElementData{"Fe", "Iron"         ,  55.845      };
      mPeriodicTable[Co] = ElementData{"Co", "Cobalt"       ,  58.933194   };
      mPeriodicTable[Ni] = ElementData{"Ni", "Nickel"       ,  58.6934     };
      mPeriodicTable[Cu] = ElementData{"Cu", "Copper"       ,  63.546      };
      mPeriodicTable[Zn] = ElementData{"Zn", "Zinc"         ,  65.38       };
      mPeriodicTable[Ga] = ElementData{"Ga", "Gallium"      ,  69.723      };
      mPeriodicTable[Ge] = ElementData{"Ge", "Germanium"    ,  72.630      };
      mPeriodicTable[As] = ElementData{"As", "Arsenic"      ,  74.921595   };
      mPeriodicTable[Se] = ElementData{"Se", "Selenium"     ,  78.971      };
      mPeriodicTable[Br] = ElementData{"Br", "Bromine"      ,  79.904      };
      mPeriodicTable[Kr] = ElementData{"Kr", "Krypton"      ,  83.798      };
      mPeriodicTable[Rb] = ElementData{"Rb", "Rubidium"     ,  85.4678     };
      mPeriodicTable[Sr] = ElementData{"Sr", "Strontium"    ,  87.62       };
      mPeriodicTable[Y ] = ElementData{"Y" , "Yttrium"      ,  88.90584    };
      mPeriodicTable[Zr] = ElementData{"Zr", "Zirconium"    ,  91.224      };
      mPeriodicTable[Nb] = ElementData{"Nb", "Niobium"      ,  92.90637    };
      mPeriodicTable[Mo] = ElementData{"Mo", "Molybdenum"   ,  95.95       };
      mPeriodicTable[Tc] = ElementData{"Tc", "Technetium"   ,  98          };
      mPeriodicTable[Ru] = ElementData{"Ru", "Ruthenium"    , 101.07       };
      mPeriodicTable[Rh] = ElementData{"Rh", "Rhodium"      , 102.90550    };
      mPeriodicTable[Pd] = ElementData{"Pd", "Palladium"    , 106.42       };
      mPeriodicTable[Ag] = ElementData{"Ag", "Silver"       , 107.8682     };
      mPeriodicTable[Cd] = ElementData{"Cd", "Cadmium"      , 112.414      };
      mPeriodicTable[In] = ElementData{"In", "Indium"       , 114.818      };
      mPeriodicTable[Sn] = ElementData{"Sn", "Tin"          , 118.710      };
      mPeriodicTable[Sb] = ElementData{"Sb", "Antimony"     , 121.760      };
      mPeriodicTable[Te] = ElementData{"Te", "Tellurium"    , 127.60       };
      mPeriodicTable[I ] = ElementData{"I" , "Iodine"       , 126.90447    };
      mPeriodicTable[Xe] = ElementData{"Xe", "Xenon"        , 131.293      };
      mPeriodicTable[Cs] = ElementData{"Dc", "Caesium"      , 132.90545196 };
      mPeriodicTable[Ba] = ElementData{"Ba", "Barium"       , 137.327      };
      mPeriodicTable[La] = ElementData{"La", "Lanthanum"    , 138.90547    };
      mPeriodicTable[Ce] = ElementData{"Ce", "Cerium"       , 140.116      };
      mPeriodicTable[Pr] = ElementData{"Pr", "Praseodymium" , 140.90766    };
      mPeriodicTable[Nd] = ElementData{"Nd", "Neodymium"    , 144.242      };
      mPeriodicTable[Pm] = ElementData{"Pm", "Promethium"   , 145          };
      mPeriodicTable[Sm] = ElementData{"Sm", "Samarium"     , 150.36       };
      mPeriodicTable[Eu] = ElementData{"Eu", "Europium"     , 151.964      };
      mPeriodicTable[Gd] = ElementData{"Gd", "Gadolinium"   , 157.25       };
      mPeriodicTable[Tb] = ElementData{"Tb", "Terbium"      , 158.92535    };
      mPeriodicTable[Dy] = ElementData{"Dy", "Dysprosium"   , 162.500      };
      mPeriodicTable[Ho] = ElementData{"Ho", "Holmium"      , 164.93033    };
      mPeriodicTable[Er] = ElementData{"Er", "Erbium"       , 167.259      };
      mPeriodicTable[Tm] = ElementData{"Tm", "Thulium"      , 168.93422    };
      mPeriodicTable[Yb] = ElementData{"Yb", "Ytterbium"    , 173.045      };
      mPeriodicTable[Lu] = ElementData{"Lu", "Lutetium"     , 174.9668     };
      mPeriodicTable[Hf] = ElementData{"Hf", "Hafnium"      , 178.49       };
      mPeriodicTable[Ta] = ElementData{"Ta", "Tantalum"     , 180.94788    };
      mPeriodicTable[W ] = ElementData{"W" , "Tungsten"     , 183.84       };
      mPeriodicTable[Re] = ElementData{"Re", "Rhenium"      , 186.207      };
      mPeriodicTable[Os] = ElementData{"Os", "Osmium"       , 190.23       };
      mPeriodicTable[Ir] = ElementData{"Ir", "Iridium"      , 192.217      };
      mPeriodicTable[Pt] = ElementData{"Pt", "Platinum"     , 195.084      };
      mPeriodicTable[Au] = ElementData{"Au", "Gold"         , 196.966569   };
      mPeriodicTable[Hg] = ElementData{"Hg", "Mercury"      , 200.592      };
      mPeriodicTable[Tl] = ElementData{"Tl", "Thallium"     , 204.38       };
      mPeriodicTable[Pb] = ElementData{"Pb", "Lead"         , 207.2        };
      mPeriodicTable[Bi] = ElementData{"Bi", "Bismuth"      , 208.98040    };
      mPeriodicTable[Po] = ElementData{"Po", "Polonium"     , 209          };
      mPeriodicTable[At] = ElementData{"At", "Astatine"     , 210          };
      mPeriodicTable[Rn] = ElementData{"Rn", "Radon"        , 222          };
      mPeriodicTable[Fr] = ElementData{"Fr", "Francium"     , 223          };
      mPeriodicTable[Ra] = ElementData{"Ra", "Radium"       , 226          };
      mPeriodicTable[Ac] = ElementData{"Ac", "Actinium"     , 227          };
      mPeriodicTable[Th] = ElementData{"Th", "Thorium"      , 232.0377     };
      mPeriodicTable[Pa] = ElementData{"Pa", "Protactinium" , 231.03588    };
      mPeriodicTable[U ] = ElementData{"U" , "Uranium"      , 238.02891    };
      mPeriodicTable[Np] = ElementData{"Np", "Neptunium"    , 237          };
      mPeriodicTable[Pu] = ElementData{"Pu", "Plutonium"    , 244          };
      mPeriodicTable[Am] = ElementData{"Am", "Americium"    , 243          };
      mPeriodicTable[Cm] = ElementData{"Cm", "Curium"       , 247          };
      mPeriodicTable[Bk] = ElementData{"Bk", "Berkelium"    , 247          };
      mPeriodicTable[Cf] = ElementData{"Cf", "Californium"  , 251          };
      mPeriodicTable[Es] = ElementData{"Es", "Einsteinium"  , 252          };
      mPeriodicTable[Fm] = ElementData{"Fm", "Fermium"      , 257          };
      mPeriodicTable[Md] = ElementData{"Md", "Mendelevium"  , 258          };
      mPeriodicTable[No] = ElementData{"No", "Nobelium"     , 259          };
      mPeriodicTable[Lr] = ElementData{"Lr", "Lawrencium"   , 266          };
      mPeriodicTable[Rf] = ElementData{"Rf", "Rutherfordium", 267          };
      mPeriodicTable[Db] = ElementData{"Db", "Dubnium"      , 268          };
      mPeriodicTable[Sg] = ElementData{"Sg", "Seaborgium"   , 269          };
      mPeriodicTable[Bh] = ElementData{"Bh", "Bohrium"      , 270          };
      mPeriodicTable[Hs] = ElementData{"Hs", "Hassium"      , 277          };
      mPeriodicTable[Mt] = ElementData{"Mt", "Meitnerium"   , 278          };
      mPeriodicTable[Ds] = ElementData{"Ds", "Darmstadtium" , 281          };
      mPeriodicTable[Rg] = ElementData{"Rg", "Roentgenium"  , 282          };
      mPeriodicTable[Cn] = ElementData{"Cn", "Copernicium"  , 285          };
      mPeriodicTable[Nh] = ElementData{"Nh", "Nihonium"     , 286          };
      mPeriodicTable[Fl] = ElementData{"Fl", "Flerovium"    , 289          };
      mPeriodicTable[Mc] = ElementData{"Mc", "Moscovium"    , 290          };
      mPeriodicTable[Lv] = ElementData{"Lv", "Livermorium"  , 293          };
      mPeriodicTable[Ts] = ElementData{"Ts", "Tennessine"   , 294          };
      mPeriodicTable[Og] = ElementData{"Og", "Oganesson"    , 294          };

}
