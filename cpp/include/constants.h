#ifndef NUMERICAL_AND_MATHEMATICAL_CONSTANTS_H
#define NUMERICAL_AND_MATHEMATICAL_CONSTANTS_H

/// Provide numerical constants.
namespace constants {
template<class Number=double>
inline constexpr Number pi()
{
    return 3.14159265358979323846264338327950288419716939937510;
}

template<class Number=double>
inline constexpr Number pion_mass()
    /// Return the pion mass in the isospin limit in MeV.
{
    return 139.57;
}

template<class Number=double>
inline constexpr Number pion_mass_squared()
    /// @brief Return the squared pion mass in the isospin limit
    /// in \f$\text{MeV}^2\f$.
{
    return pion_mass<Number>()*pion_mass<Number>();
}

/// PDG (Particle Data Group) values
namespace pdg {
constexpr double pion_decay{92.28}; // in MeV
constexpr double charged_pion_mass{139.57061}; // in MeV
constexpr double omega_mass{782.65}; // in MeV
} // pdg

/// @brief FLAG (Flavour Lattice Averaging Group) values (2019 edition)
namespace flag {
/// \f$F_\pi/F\f$ for \f$N_f = 2+1+1\f$.
constexpr double pion_decay_ratio_four_flavours{1.077};

/// \f$F_\pi/F\f$ for \f$N_f = 2+1\f$.
constexpr double pion_decay_ratio_three_flavours{1.062};

/// \f$F_\pi/F\f$ for \f$N_f = 2\f$.
constexpr double pion_decay_ratio_two_flavours{1.073};
} // flag

constexpr double pion_decay_chiral{
    pdg::pion_decay/flag::pion_decay_ratio_three_flavours};
} // constants


#endif // NUMERICAL_AND_MATHEMATICAL_CONSTANTS_H
