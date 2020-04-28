#ifndef OMNES_CURVE_h
#define OMNES_CURVE_h

#include "grid.h"
#include "omnes.h"
#include "type_aliases.h"

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

/// Omnes function with cut along a somewhat general curve.
namespace curved_omnes {
using type_aliases::Complex;
using type_aliases::CFunction;

std::vector<Complex> first_points(const grid::Curve& curve, std::size_t size);
    ///< Extract the first `size` boundary points of `curve`.

bool on_second_sheet(const std::vector<Complex>& points,
        const Complex& mandelstam_s);
    ///< Determine if `mandelstam_s` is on the second sheet.

/// An Omnes function with a cut along `curve`.
class CurvedOmnes {
public:
    template<typename C>
    CurvedOmnes(omnes::OmnesF o, CFunction amplitude, C&& curve)
        /// @param o The omnes function with the usual right-hand cut.
        /// @param amplitude The two-to-two particle scattering amplitude
        /// associated with the phase of `o`.
        /// @param curve The branch cut of the returned Omnes function. Currently,
        /// this needs to be a piecewise curve with at least four knots, where
        /// the first four knots form a rectangle extending into the lower half
        /// plane. This is due to the non-general inner-workings of
        /// `on_second_sheet`.
        : o{std::move(o)}, amplitude{std::move(amplitude)},
        points{first_points(std::forward<C>(curve), 4)}
    {
        static_assert(std::is_base_of_v<grid::Curve,
                                        std::remove_reference_t<C>>,
            "Cut needs to be specified as a class inheriting from `Curve`.");
    }

    Complex operator()(const Complex& mandelstam_s) const
    {
        if (on_second_sheet(points, mandelstam_s))
            return omnes::second_sheet(o, amplitude, mandelstam_s);
        return o(mandelstam_s);
    }
private:
    omnes::OmnesF o;
    CFunction amplitude;
    std::vector<Complex> points;
};
} // curved_omnes

#endif // OMNES_CURVE_h
