#include "curved_omnes.h"

namespace curved_omnes {
std::vector<Complex> first_points(const grid::Curve& curve, std::size_t size)
{
    const std::vector<double> boundaries{curve.boundaries()};
    if (boundaries.size() < size) {
        std::stringstream message;
        message << "Tried to retrieve " << size << " elements, but curve has "
            << "only " << boundaries.size() << " boundary points.";
        throw std::domain_error{message.str()};
    }
    std::vector<Complex> points;
    points.reserve(size);
    const auto start{boundaries.cbegin()};
    std::transform(start, start + size, std::back_inserter(points),
            [&curve](double x){return curve.curve_func(x);});
    return points;
}

bool on_second_sheet(const std::vector<Complex>& points,
        const Complex& mandelstam_s)
{
    return points[0].real() < mandelstam_s.real()
        && mandelstam_s.real() < points[3].real()
        && points[1].imag() < mandelstam_s.imag()
        && mandelstam_s.imag() < 0.0;
}
} // curved_omnes
