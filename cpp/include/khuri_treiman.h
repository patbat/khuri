#ifndef KHURI_TREIMAN_H
#define KHURI_TREIMAN_H

#include "kernel.h"
#include "piecewise.h"

/// @brief Solve Khuri-Treiman (KT) equations for the scattering/decay involving
/// three pions with arbitrary mass and one particle with I=0, J=1, P=C=-1
/// and arbitrary mass.
///
/// The equations are solved via the modified Gasser-Rusetsky method.
namespace khuri_treiman {
using grid::Curve;
using grid::Point;
using grid::Grid;

using piecewise::Piecewise;
using piecewise::Adaptive;
using piecewise::Vector_decay;
using piecewise::Real;

using kernel::Basis;
using kernel::Complex;
using kernel::CFunction;
using kernel::Method;
} // khuri_treiman

#endif // KHURI_TREIMAN_H
