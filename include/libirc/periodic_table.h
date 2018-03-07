#ifndef IRC_PERIODIC_TABLE_H
#define IRC_PERIODIC_TABLE_H

#include "conversion.h"

#include <array>
#include <iostream>
#include <string>

namespace irc {

namespace periodic_table {

/// Number of atomic species in the periodic table
constexpr size_t pt_size{96};

/// Atomic symbols
constexpr std::array<char[3], pt_size> symbols = {
    {"",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
     "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
     "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
     "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
     "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
     "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
     "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
     "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am"}};

/// Atomic masses
constexpr std::array<double, pt_size> masses = {{
    0.0000,   //   Ghost Atom      X       0
    1.0079,   //   Hydrogen        H       1
    4.0026,   //   Helium          He      2
    6.941,    //   Lithium         Li      3
    9.0122,   //   Beryllium       Be      4
    10.811,   //   Boron           B       5
    12.0107,  //   Carbon          C       6
    14.0067,  //   Nitrogen        N       7
    15.9994,  //   Oxygen          O       8
    18.9984,  //   Fluorine        F       9
    20.1797,  //   Neon            Ne      10
    22.9897,  //   Sodium          Na      11
    24.305,   //   Magnesium       Mg      12
    26.9815,  //   Aluminum        Al      13
    28.0855,  //   Silicon         Si      14
    30.9738,  //   Phosphorus      P       15
    32.065,   //   Sulfur          S       16
    35.453,   //   Chlorine        Cl      17
    39.948,   //   Argon           Ar      18
    39.0983,  //   Potassium       K       19
    40.078,   //   Calcium         Ca      20
    44.9559,  //   Scandium        Sc      21
    47.867,   //   Titanium        Ti      22
    50.9415,  //   Vanadium        V       23
    51.9961,  //   Chromium        Cr      24
    54.938,   //   Manganese       Mn      25
    55.845,   //   Iron            Fe      26
    58.9332,  //   Cobalt          Co      27
    58.6934,  //   Nickel          Ni      28
    63.546,   //   Copper          Cu      29
    65.39,    //   Zinc            Zn      30
    69.723,   //   Gallium         Ga      31
    72.64,    //   Germanium       Ge      32
    74.9216,  //   Arsenic         As      33
    78.96,    //   Selenium        Se      34
    79.904,   //   Bromine         Br      35
    83.8,     //   Krypton         Kr      36
    85.4678,  //   Rubidium        Rb      37
    87.62,    //   Strontium       Sr      38
    88.9059,  //   Yttrium         Y       39
    91.224,   //   Zirconium       Zr      40
    92.9064,  //   Niobium         Nb      41
    95.94,    //   Molybdenum      Mo      42
    98.,      //   Technetium      Tc      43
    101.07,   //   Ruthenium       Ru      44
    102.9055, //   Rhodium         Rh      45
    106.42,   //   Palladium       Pd      46
    107.8682, //   Silver          Ag      47
    112.411,  //   Cadmium         Cd      48
    114.818,  //   Indium          In      49
    118.71,   //   Tin             Sn      50
    121.76,   //   Antimony        Sb      51
    127.6,    //   Tellurium       Te      52
    126.9045, //   Iodine          I       53
    131.293,  //   Xenon           Xe      54
    132.9055, //   Cesium          Cs      55
    137.327,  //   Barium          Ba      56
    138.9055, //   Lanthanum       La      57
    140.116,  //   Cerium          Ce      58
    140.9077, //   Praseodymium    Pr      59
    144.24,   //   Neodymium       Nd      60
    145,      //   Promethium      Pm      61
    150.36,   //   Samarium        Sm      62
    151.964,  //   Europium        Eu      63
    157.25,   //   Gadolinium      Gd      64
    158.9253, //   Terbium         Tb      65
    162.5,    //   Dysprosium      Dy      66
    164.9303, //   Holmium         Ho      67
    167.259,  //   Erbium          Er      68
    168.9342, //   Thulium         Tm      69
    173.04,   //   Ytterbium       Yb      70
    174.967,  //   Lutetium        Lu      71
    178.49,   //   Hafnium         Hf      72
    180.9479, //   Tantalum        Ta      73
    183.84,   //   Tungsten        W       74
    186.207,  //   Rhenium         Re      75
    190.23,   //   Osmium          Os      76
    192.217,  //   Iridium         Ir      77
    195.078,  //   Platinum        Pt      78
    196.9665, //   Gold            Au      79
    200.59,   //   Mercury         Hg      80
    204.3833, //   Thallium        Tl      81
    207.2,    //   Lead            Pb      82
    208.9804, //   Bismuth         Bi      83
    209.,     //   Polonium        Po      84
    210.,     //   Astatine        At      85
    222.,     //   Radon           Rn      86
    223.,     //   Francium        Fr      87
    226.,     //   Radium          Ra      88
    227.,     //   Actinium        Ac      89
    232.0381, //   Thorium         Th      90
    231.0359, //   Protactinium    Pa      91
    238.0289, //   Uranium         U       92
    237.,     //   Neptunium       Np      93
    244.,     //   Plutonium       Pu      94
    243.,     //   Americium       Am      95
}};

using tools::conversion::angstrom_to_bohr;

/// Covalent radii
///
/// J. C. Slater, "Atomic Radii in Crystals",
/// J. Phys. Chem. 41, 3199-3204 (1964)
constexpr std::array<double, pt_size> covalent_radii = {{
    0.00 * angstrom_to_bohr, //   Ghost Atom      X       0
    0.25 * angstrom_to_bohr, //   Hydrogen        H       1
    0.00 * angstrom_to_bohr, //   Helium          He      2
    1.45 * angstrom_to_bohr, //   Lithium         Li      3
    1.05 * angstrom_to_bohr, //   Beryllium       Be      4
    0.85 * angstrom_to_bohr, //   Boron           B       5
    0.70 * angstrom_to_bohr, //   Carbon          C       6
    0.65 * angstrom_to_bohr, //   Nitrogen        N       7
    0.60 * angstrom_to_bohr, //   Oxygen          O       8
    0.50 * angstrom_to_bohr, //   Fluorine        F       9
    0.00 * angstrom_to_bohr, //   Neon            Ne      10
    1.80 * angstrom_to_bohr, //   Sodium          Na      11
    1.50 * angstrom_to_bohr, //   Magnesium       Mg      12
    1.25 * angstrom_to_bohr, //   Aluminum        Al      13
    1.10 * angstrom_to_bohr, //   Silicon         Si      14
    1.00 * angstrom_to_bohr, //   Phosphorus      P       15
    1.00 * angstrom_to_bohr, //   Sulfur          S       16
    1.00 * angstrom_to_bohr, //   Chlorine        Cl      17
    0.00 * angstrom_to_bohr, //   Argon           Ar      18
    2.20 * angstrom_to_bohr, //   Potassium       K       19
    1.80 * angstrom_to_bohr, //   Calcium         Ca      20
    1.60 * angstrom_to_bohr, //   Scandium        Sc      21
    1.40 * angstrom_to_bohr, //   Titanium        Ti      22
    1.35 * angstrom_to_bohr, //   Vanadium        V       23
    1.40 * angstrom_to_bohr, //   Chromium        Cr      24
    1.40 * angstrom_to_bohr, //   Manganese       Mn      25
    1.40 * angstrom_to_bohr, //   Iron            Fe      26
    1.35 * angstrom_to_bohr, //   Cobalt          Co      27
    1.35 * angstrom_to_bohr, //   Nickel          Ni      28
    1.35 * angstrom_to_bohr, //   Copper          Cu      29
    1.35 * angstrom_to_bohr, //   Zinc            Zn      30
    1.30 * angstrom_to_bohr, //   Gallium         Ga      31
    1.25 * angstrom_to_bohr, //   Germanium       Ge      32
    1.15 * angstrom_to_bohr, //   Arsenic         As      33
    1.15 * angstrom_to_bohr, //   Selenium        Se      34
    1.15 * angstrom_to_bohr, //   Bromine         Br      35
    0.00 * angstrom_to_bohr, //   Krypton         Kr      36
    2.35 * angstrom_to_bohr, //   Rubidium        Rb      37
    2.00 * angstrom_to_bohr, //   Strontium       Sr      38
    1.80 * angstrom_to_bohr, //   Yttrium         Y       39
    1.55 * angstrom_to_bohr, //   Zirconium       Zr      40
    1.45 * angstrom_to_bohr, //   Niobium         Nb      41
    1.45 * angstrom_to_bohr, //   Molybdenum      Mo      42
    1.35 * angstrom_to_bohr, //   Technetium      Tc      43
    1.30 * angstrom_to_bohr, //   Ruthenium       Ru      44
    1.35 * angstrom_to_bohr, //   Rhodium         Rh      45
    1.40 * angstrom_to_bohr, //   Palladium       Pd      46
    1.60 * angstrom_to_bohr, //   Silver          Ag      47
    1.55 * angstrom_to_bohr, //   Cadmium         Cd      48
    1.55 * angstrom_to_bohr, //   Indium          In      49
    1.45 * angstrom_to_bohr, //   Tin             Sn      50
    1.45 * angstrom_to_bohr, //   Antimony        Sb      51
    1.40 * angstrom_to_bohr, //   Tellurium       Te      52
    1.40 * angstrom_to_bohr, //   Iodine          I       53
    0.00 * angstrom_to_bohr, //   Xenon           Xe      54
    2.60 * angstrom_to_bohr, //   Cesium          Cs      55
    2.15 * angstrom_to_bohr, //   Barium          Ba      56
    1.95 * angstrom_to_bohr, //   Lanthanum       La      57
    1.85 * angstrom_to_bohr, //   Cerium          Ce      58
    1.85 * angstrom_to_bohr, //   Praseodymium    Pr      59
    1.85 * angstrom_to_bohr, //   Neodymium       Nd      60
    1.85 * angstrom_to_bohr, //   Promethium      Pm      61
    1.85 * angstrom_to_bohr, //   Samarium        Sm      62
    1.85 * angstrom_to_bohr, //   Europium        Eu      63
    1.80 * angstrom_to_bohr, //   Gadolinium      Gd      64
    1.75 * angstrom_to_bohr, //   Terbium         Tb      65
    1.75 * angstrom_to_bohr, //   Dysprosium      Dy      66
    1.75 * angstrom_to_bohr, //   Holmium         Ho      67
    1.75 * angstrom_to_bohr, //   Erbium          Er      68
    1.75 * angstrom_to_bohr, //   Thulium         Tm      69
    1.75 * angstrom_to_bohr, //   Ytterbium       Yb      70
    1.75 * angstrom_to_bohr, //   Lutetium        Lu      71
    1.55 * angstrom_to_bohr, //   Hafnium         Hf      72
    1.45 * angstrom_to_bohr, //   Tantalum        Ta      73
    1.35 * angstrom_to_bohr, //   Tungsten        W       74
    1.35 * angstrom_to_bohr, //   Rhenium         Re      75
    1.30 * angstrom_to_bohr, //   Osmium          Os      76
    1.35 * angstrom_to_bohr, //   Iridium         Ir      77
    1.35 * angstrom_to_bohr, //   Platinum        Pt      78
    1.35 * angstrom_to_bohr, //   Gold            Au      79
    1.50 * angstrom_to_bohr, //   Mercury         Hg      80
    1.90 * angstrom_to_bohr, //   Thallium        Tl      81
    0.00 * angstrom_to_bohr, //   Lead            Pb      82
    1.60 * angstrom_to_bohr, //   Bismuth         Bi      83
    1.60 * angstrom_to_bohr, //   Polonium        Po      84
    0.00,                    //   Astatine        At      85
    0.00 * angstrom_to_bohr, //   Radon           Rn      86
    0.00,                    //   Francium        Fr      87
    2.15 * angstrom_to_bohr, //   Radium          Ra      88
    1.95 * angstrom_to_bohr, //   Actinium        Ac      89
    1.80 * angstrom_to_bohr, //   Thorium         Th      90
    1.80 * angstrom_to_bohr, //   Protactinium    Pa      91
    1.75 * angstrom_to_bohr, //   Uranium         U       92
    1.75 * angstrom_to_bohr, //   Neptunium       Np      93
    1.75 * angstrom_to_bohr, //   Plutonium       Pu      94
    1.75 * angstrom_to_bohr, //   Americium       Am      95
}};

/// Van der Waals radii
///
/// S. S. Batsanov, "Van der Waals Radii of Elements",
/// Inorganic Materials 37, 871-855 (2001)
constexpr std::array<double, pt_size> vdw_radii = {{
    0.00 * angstrom_to_bohr, //   Ghost Atom      X       0
    1.20 * angstrom_to_bohr, //   Hydrogen        H       1
    0 * angstrom_to_bohr,    //   Helium          He      2
    2.24 * angstrom_to_bohr, //   Lithium         Li      3
    1.86 * angstrom_to_bohr, //   Beryllium       Be      4
    1.74 * angstrom_to_bohr, //   Boron           B       5
    1.70 * angstrom_to_bohr, //   Carbon          C       6
    1.55 * angstrom_to_bohr, //   Nitrogen        N       7
    1.52 * angstrom_to_bohr, //   Oxygen          O       8
    1.47 * angstrom_to_bohr, //   Fluorine        F       9
    0 * angstrom_to_bohr,    //   Neon            Ne      10
    2.57 * angstrom_to_bohr, //   Sodium          Na      11
    2.27 * angstrom_to_bohr, //   Magnesium       Mg      12
    2.11 * angstrom_to_bohr, //   Aluminum        Al      13
    2.06 * angstrom_to_bohr, //   Silicon         Si      14
    0 * angstrom_to_bohr,    //   Phosphorus      P       15
    1.80 * angstrom_to_bohr, //   Sulfur          S       16
    1.75 * angstrom_to_bohr, //   Chlorine        Cl      17
    0 * angstrom_to_bohr,    //   Argon           Ar      18
    3.00 * angstrom_to_bohr, //   Potassium       K       19
    2.61 * angstrom_to_bohr, //   Calcium         Ca      20
    2.28 * angstrom_to_bohr, //   Scandium        Sc      21
    2.14 * angstrom_to_bohr, //   Titanium        Ti      22
    2.03 * angstrom_to_bohr, //   Vanadium        V       23
    1.97 * angstrom_to_bohr, //   Chromium        Cr      24
    1.96 * angstrom_to_bohr, //   Manganese       Mn      25
    1.96 * angstrom_to_bohr, //   Iron            Fe      26
    1.95 * angstrom_to_bohr, //   Cobalt          Co      27
    1.94 * angstrom_to_bohr, //   Nickel          Ni      28
    2.00 * angstrom_to_bohr, //   Copper          Cu      29
    2.02 * angstrom_to_bohr, //   Zinc            Zn      30
    2.08 * angstrom_to_bohr, //   Gallium         Ga      31
    2.13 * angstrom_to_bohr, //   Germanium       Ge      32
    2.16 * angstrom_to_bohr, //   Arsenic         As      33
    0 * angstrom_to_bohr,    //   Selenium        Se      34
    1.85 * angstrom_to_bohr, //   Bromine         Br      35
    0 * angstrom_to_bohr,    //   Krypton         Kr      36
    0 * angstrom_to_bohr,    //   Rubidium        Rb      37
    2.78 * angstrom_to_bohr, //   Strontium       Sr      38
    2.45 * angstrom_to_bohr, //   Yttrium         Y       39
    2.25 * angstrom_to_bohr, //   Zirconium       Zr      40
    2.13 * angstrom_to_bohr, //   Niobium         Nb      41
    2.06 * angstrom_to_bohr, //   Molybdenum      Mo      42
    2.04 * angstrom_to_bohr, //   Technetium      Tc      43
    2.02 * angstrom_to_bohr, //   Ruthenium       Ru      44
    2.02 * angstrom_to_bohr, //   Rhodium         Rh      45
    2.05 * angstrom_to_bohr, //   Palladium       Pd      46
    2.13 * angstrom_to_bohr, //   Silver          Ag      47
    2.17 * angstrom_to_bohr, //   Cadmium         Cd      48
    2.24 * angstrom_to_bohr, //   Indium          In      49
    2.29 * angstrom_to_bohr, //   Tin             Sn      50
    2.33 * angstrom_to_bohr, //   Antimony        Sb      51
    0 * angstrom_to_bohr,    //   Tellurium       Te      52
    1.98 * angstrom_to_bohr, //   Iodine          I       53
    0 * angstrom_to_bohr,    //   Xenon           Xe      54
    0 * angstrom_to_bohr,    //   Cesium          Cs      55
    2.85 * angstrom_to_bohr, //   Barium          Ba      56
    2.51 * angstrom_to_bohr, //   Lanthanum       La      57
    0 * angstrom_to_bohr,    //   Cerium          Ce      58
    0 * angstrom_to_bohr,    //   Praseodymium    Pr      59
    0 * angstrom_to_bohr,    //   Neodymium       Nd      60
    0 * angstrom_to_bohr,    //   Promethium      Pm      61
    0 * angstrom_to_bohr,    //   Samarium        Sm      62
    0 * angstrom_to_bohr,    //   Europium        Eu      63
    0 * angstrom_to_bohr,    //   Gadolinium      Gd      64
    0 * angstrom_to_bohr,    //   Terbium         Tb      65
    0 * angstrom_to_bohr,    //   Dysprosium      Dy      66
    0 * angstrom_to_bohr,    //   Holmium         Ho      67
    0 * angstrom_to_bohr,    //   Erbium          Er      68
    0 * angstrom_to_bohr,    //   Thulium         Tm      69
    0 * angstrom_to_bohr,    //   Ytterbium       Yb      70
    0 * angstrom_to_bohr,    //   Lutetium        Lu      71
    2.24 * angstrom_to_bohr, //   Hafnium         Hf      72
    0 * angstrom_to_bohr,    //   Tantalum        Ta      73
    0 * angstrom_to_bohr,    //   Tungsten        W       74
    0 * angstrom_to_bohr,    //   Rhenium         Re      75
    2.03 * angstrom_to_bohr, //   Osmium          Os      76
    2.03 * angstrom_to_bohr, //   Iridium         Ir      77
    2.06 * angstrom_to_bohr, //   Platinum        Pt      78
    2.13 * angstrom_to_bohr, //   Gold            Au      79
    2.17 * angstrom_to_bohr, //   Mercury         Hg      80
    0 * angstrom_to_bohr,    //   Thallium        Tl      81
    2.36 * angstrom_to_bohr, //   Lead            Pb      82
    0 * angstrom_to_bohr,    //   Bismuth         Bi      83
    0 * angstrom_to_bohr,    //   Polonium        Po      84
    0 * angstrom_to_bohr,    //   Astatine        At      85
    0 * angstrom_to_bohr,    //   Radon           Rn      86
    0 * angstrom_to_bohr,    //   Francium        Fr      87
    0 * angstrom_to_bohr,    //   Radium          Ra      88
    0 * angstrom_to_bohr,    //   Actinium        Ac      89
    2.43 * angstrom_to_bohr, //   Thorium         Th      90
    0 * angstrom_to_bohr,    //   Protactinium    Pa      91
    2.17 * angstrom_to_bohr, //   Uranium         U       92
    0 * angstrom_to_bohr,    //   Neptunium       Np      93
    0 * angstrom_to_bohr,    //   Plutonium       Pu      94
    0 * angstrom_to_bohr,    //   Americium       Am      95
}};

/// Check if atomic number \param an is valid
///
/// \param an Atomic number
bool valid_atomic_number(size_t an) noexcept {
  return an > 0 && an < pt_size;
}

/// Return atomic number corresponding to symbol
///
/// \param symbol Atomic symbol
/// \return Atomic number
size_t atomic_number(const std::string& symbol) {
  size_t an{0};
  
  for (size_t i = 0; i < pt_size; i++) {
    if (symbol == symbols[i]) {
      an = i;
      break;
    }
  }
  
  if (an == 0) {
    throw std::logic_error("Invalid atomic symbol.");
  }
  
  return an;
}

} // namespace periodic_table

} // namespace irc

#endif // IRC_PERIODIC_TABLE_H
