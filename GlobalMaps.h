#ifndef _GLOBAL_MAPS_TYPE__
#define _GLOBAL_MAPS_TYPE__

#include <map>
#include <string>

using namespace std;

map< string, ElementType > gElementTypeMap = {
  {"e", e},  //  electron  0.0
  {"H", H},  // 	Hydrogen 1.008
  {"He", He}, // 	Helium 4.002602
  {"Li", Li}, // 	Lithium 6.94
  {"Be", Be}, // 	Beryllium 9.0121831
  {"B", B},  // 	Boron 10.81
  {"C", C},  // 	Carbon 12.011
  {"N", N},  // 	Nitrogen 14.007
  {"O", O},  // 	Oxygen 15.999
  {"F", F},  // 	Fluorine 18.998403163
  {"Ne", Ne}, // 	Neon 20.1797
  {"Na", Na}, // 	Sodium 22.98976928
  {"Mg", Mg}, // 	Magnesium 24.305
  {"Al", Al}, // 	Aluminium 26.9815385
  {"Si", Si}, // 	Silicon 28.085
  {"P", P},  // 	Phosphorus 30.973761998
  {"S", S},  // 	Sulfur 32.06
  {"Cl", Cl}, // 	Chlorine 35.45
  {"Ar", Ar}, // 	Argon 39.948
  {"K", K},  // 	Potassium 39.0983
  {"Ca", Ca}, // 	Calcium 40.078
  {"Sc", Sc}, // 	Scandium 44.955908
  {"Ti", Ti}, // 	Titanium 47.867
  {"V", V},  // 	Vanadium 50.9415
  {"Cr", Cr}, // 	Chromium 51.9961
  {"Mn", Mn}, // 	Manganese 54.938044
  {"Fe", Fe}, // 	Iron 55.845
  {"Co", Co}, // 	Cobalt 58.933194
  {"Ni", Ni}, // 	Nickel 58.6934
  {"Cu", Cu}, // 	Copper 63.546
  {"Zn", Zn}, // 	Zinc 65.38
  {"Ga", Ga}, // 	Gallium 69.723
  {"Ge", Ge}, // 	Germanium 72.630
  {"As", As}, // 	Arsenic 74.921595
  {"Se", Se}, // 	Selenium 78.971
  {"Br", Br}, // 	Bromine 79.904
  {"Kr", Kr}, // 	Krypton 83.798
  {"Rb", Rb}, // 	Rubidium 85.4678
  {"Sr", Sr}, // 	Strontium 87.62
  {"Y", Y},  // 	Yttrium 88.90584
  {"Zr", Zr}, // 	Zirconium 91.224
  {"Nb", Nb}, // 	Niobium 92.90637
  {"Mo", Mo}, // 	Molybdenum 95.95
  {"Tc", Tc}, // 	Technetium 98
  {"Ru", Ru}, // 	Ruthenium 101.07
  {"Rh", Rh}, // 	Rhodium 102.90550
  {"Pd", Pd}, // 	Palladium 106.42
  {"Ag", Ag}, // 	Silver 107.8682
  {"Cd", Cd}, //	Cadmium 112.414
  {"In", In}, // 	Indium 114.818
  {"Sn", Sn}, // 	Tin 118.710
  {"Sb", Sb}, // 	Antimony 121.760
  {"Te", Te}, // 	Tellurium 127.60
  {"I", I},  // 	Iodine 126.90447
  {"Xe", Xe}, // 	Xenon 131.293
  {"Cs", Cs}, // 	Caesium 132.90545196
  {"Ba", Ba}, // 	Barium 137.327
  {"La", La}, // 	Lanthanum 138.90547
  {"Ce", Ce}, // 	Cerium 140.116
  {"Pr", Pr}, // 	Praseodymium 140.90766
  {"Nd", Nd}, // 	Neodymium 144.242
  {"Pm", Pm}, // 	Promethium 145
  {"Sm", Sm}, // 	Samarium 150.36
  {"Eu", Eu}, // 	Europium 151.964
  {"Gd", Gd}, // 	Gadolinium 157.25
  {"Tb", Tb}, // 	Terbium 158.92535
  {"Dy", Dy}, // 	Dysprosium 162.500
  {"Ho", Ho}, // 	Holmium 164.93033
  {"Er", Er}, // 	Erbium 167.259
  {"Tm", Tm}, // 	Thulium 168.93422
  {"Yb", Yb}, // 	Ytterbium 173.045
  {"Lu", Lu}, // 	Lutetium 174.9668
  {"Hf", Hf}, // 	Hafnium 178.49
  {"Ta", Ta}, // 	Tantalum 180.94788
  {"W", W},  // 	Tungsten 183.84
  {"Re", Re}, // 	Rhenium 186.207
  {"Os", Os}, // 	Osmium 190.23
  {"Ir", Ir}, // 	Iridium 192.217
  {"Pt", Pt}, // 	Platinum 195.084
  {"Au", Au}, // 	Gold 196.966569
  {"Hg", Hg}, // 	Mercury 200.592
  {"Tl", Tl}, // 	Thallium 204.38
  {"Pb", Pb}, // 	Lead 207.2
  {"Bi", Bi}, // 	Bismuth 208.98040
  {"Po", Po}, // 	Polonium 209
  {"At", At}, // 	Astatine 210
  {"Rn", Rn}, // 	Radon 222
  {"Fr", Fr}, // 	Francium 223
  {"Ra", Ra}, // 	Radium 226
  {"Ac", Ac}, // 	Actinium 227
  {"Th", Th}, // 	Thorium 232.0377
  {"Pa", Pa}, // 	Protactinium 231.03588
  {"U", U},  // 	Uranium 238.02891
  {"Np", Np}, // 	Neptunium 237
  {"Pu", Pu}, // 	Plutonium 244
  {"Am", Am}, // 	Americium 243
  {"Cm", Cm}, // 	Curium 247
  {"Bk", Bk}, // 	Berkelium 247
  {"Cf", Cf}, // 	Californium 251
  {"Es", Es}, // 	Einsteinium 252
  {"Fm", Fm}, // 	Fermium 257
  {"Md", Md}, // 	Mendelevium 258
  {"No", No}, // 	Nobelium 259
  {"Lr", Lr}, // 	Lawrencium 266
  {"Rf", Rf}, // 	Rutherfordium 267
  {"Db", Db}, // 	Dubnium 268
  {"Sg", Sg}, // 	Seaborgium 269
  {"Bh", Bh}, // 	Bohrium 270
  {"Hs", Hs}, // 	Hassium 277
  {"Mt", Mt}, // 	Meitnerium 278
  {"Ds", Ds}, // 	Darmstadtium 281
  {"Rg", Rg}, // 	Roentgenium 282
  {"Cn", Cn}, // 	Copernicium 285
  {"Nh", Nh}, // 	Nihonium 286
  {"Fl", Fl}, // 	Flerovium 289
  {"Mc", Mc}, // 	Moscovium 290
  {"Lv", Lv}, // 	Livermorium 293
  {"Ts", Ts}, // 	Tennessine 294
  {"Og", Og}, // 	Oganesson 294
  {"Unknown", ElementSize}
};

#endif
