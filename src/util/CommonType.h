#ifndef ____CommonType__
#define ____CommonType__

#include <map>
#include <string>

using namespace std;

namespace CommonType {
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
    UNDEFINED,
    ALPHA_CARBON,
    MASS_WEIGHTED_MEAN
  };

  enum ElementType{
    e,  //  electron  0.0
    H,  // 	Hydrogen 1.008
    He, // 	Helium 4.002602
    Li, // 	Lithium 6.94
    Be, // 	Beryllium 9.0121831
    B,  // 	Boron 10.81
    C,  // 	Carbon 12.011
    N,  // 	Nitrogen 14.007
    O,  // 	Oxygen 15.999
    F,  // 	Fluorine 18.998403163
    Ne, // 	Neon 20.1797
    Na, // 	Sodium 22.98976928
    Mg, // 	Magnesium 24.305
    Al, // 	Aluminium 26.9815385
    Si, // 	Silicon 28.085
    P,  // 	Phosphorus 30.973761998
    S,  // 	Sulfur 32.06
    Cl, // 	Chlorine 35.45
    Ar, // 	Argon 39.948
    K,  // 	Potassium 39.0983
    Ca, // 	Calcium 40.078
    Sc, // 	Scandium 44.955908
    Ti, // 	Titanium 47.867
    V,  // 	Vanadium 50.9415
    Cr, // 	Chromium 51.9961
    Mn, // 	Manganese 54.938044
    Fe, // 	Iron 55.845
    Co, // 	Cobalt 58.933194
    Ni, // 	Nickel 58.6934
    Cu, // 	Copper 63.546
    Zn, // 	Zinc 65.38
    Ga, // 	Gallium 69.723
    Ge, // 	Germanium 72.630
    As, // 	Arsenic 74.921595
    Se, // 	Selenium 78.971
    Br, // 	Bromine 79.904
    Kr, // 	Krypton 83.798
    Rb, // 	Rubidium 85.4678
    Sr, // 	Strontium 87.62
    Y,  // 	Yttrium 88.90584
    Zr, // 	Zirconium 91.224
    Nb, // 	Niobium 92.90637
    Mo, // 	Molybdenum 95.95
    Tc, // 	Technetium 98
    Ru, // 	Ruthenium 101.07
    Rh, // 	Rhodium 102.90550
    Pd, // 	Palladium 106.42
    Ag, // 	Silver 107.8682
    Cd, //	Cadmium 112.414
    In, // 	Indium 114.818
    Sn, // 	Tin 118.710
    Sb, // 	Antimony 121.760
    Te, // 	Tellurium 127.60
    I,  // 	Iodine 126.90447
    Xe, // 	Xenon 131.293
    Cs, // 	Caesium 132.90545196
    Ba, // 	Barium 137.327
    La, // 	Lanthanum 138.90547
    Ce, // 	Cerium 140.116
    Pr, // 	Praseodymium 140.90766
    Nd, // 	Neodymium 144.242
    Pm, // 	Promethium 145
    Sm, // 	Samarium 150.36
    Eu, // 	Europium 151.964
    Gd, // 	Gadolinium 157.25
    Tb, // 	Terbium 158.92535
    Dy, // 	Dysprosium 162.500
    Ho, // 	Holmium 164.93033
    Er, // 	Erbium 167.259
    Tm, // 	Thulium 168.93422
    Yb, // 	Ytterbium 173.045
    Lu, // 	Lutetium 174.9668
    Hf, // 	Hafnium 178.49
    Ta, // 	Tantalum 180.94788
    W,  // 	Tungsten 183.84
    Re, // 	Rhenium 186.207
    Os, // 	Osmium 190.23
    Ir, // 	Iridium 192.217
    Pt, // 	Platinum 195.084
    Au, // 	Gold 196.966569
    Hg, // 	Mercury 200.592
    Tl, // 	Thallium 204.38
    Pb, // 	Lead 207.2
    Bi, // 	Bismuth 208.98040
    Po, // 	Polonium 209
    At, // 	Astatine 210
    Rn, // 	Radon 222
    Fr, // 	Francium 223
    Ra, // 	Radium 226
    Ac, // 	Actinium 227
    Th, // 	Thorium 232.0377
    Pa, // 	Protactinium 231.03588
    U,  // 	Uranium 238.02891
    Np, // 	Neptunium 237
    Pu, // 	Plutonium 244
    Am, // 	Americium 243
    Cm, // 	Curium 247
    Bk, // 	Berkelium 247
    Cf, // 	Californium 251
    Es, // 	Einsteinium 252
    Fm, // 	Fermium 257
    Md, // 	Mendelevium 258
    No, // 	Nobelium 259
    Lr, // 	Lawrencium 266
    Rf, // 	Rutherfordium 267
    Db, // 	Dubnium 268
    Sg, // 	Seaborgium 269
    Bh, // 	Bohrium 270
    Hs, // 	Hassium 277
    Mt, // 	Meitnerium 278
    Ds, // 	Darmstadtium 281
    Rg, // 	Roentgenium 282
    Cn, // 	Copernicium 285
    Nh, // 	Nihonium 286
    Fl, // 	Flerovium 289
    Mc, // 	Moscovium 290
    Lv, // 	Livermorium 293
    Ts, // 	Tennessine 294
    Og  // 	Oganesson 294
  };
}


#endif
