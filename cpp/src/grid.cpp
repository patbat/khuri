#include "grid.h"

namespace grid {
Knots generate_knots(double start, double end, std::size_t points)
{
    const gsl::Gauss_Legendre g{points};
    Knots path(points);
    for (std::size_t i{0}; i<points; ++i)
        path[i] = g.point(start,end,i);
    return path;
}
} // grid
