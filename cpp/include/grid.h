#ifndef KERNEL_GRID_H
#define KERNEL_GRID_H

#include "gsl_interface.h"

#include <algorithm>
#include <complex>
#include <iterator>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

/// Grid used in solution of KT equations.

/// Here and in the following, x refers usually to an integration variable in
/// the Mandelstam-s plane, while z is the cosine of the scattering angle.
/// Gauss-Legendre quadrature is employed in solving the integral equations.
/// That is, the integrands are sampled on a grid in the (x,z) plane (here x
/// might be complex), the precise location of the sampling points and the
/// associated weights are determined via the Gauss-Legendre method.
/// This header contains general facilites to compute such a grid, the most
/// import is the class `Grid`.
namespace grid {
using Complex = std::complex<double>;

/// (point,weight) pair
using Knot = std::pair<double,double>;
using Knots = std::vector<Knot>;

/// (point,weight,derivative) triples
template<typename T1, typename T2=double>
using Sampling_points = std::vector<std::tuple<T1,T2,T1>>;

Knots generate_knots(double start, double end, std::size_t points);
    ///< @brief Return (point,weight) pairs for Gauss-Legendre integration in
    ///< interval [`start`,`end`].

template<typename F1, typename F2>
auto knots_along_curve(double start, double end,
        std::size_t points, const F1& curve, const F2& derivative)
    /// Compute `curve` and `derivative` at Gauss-Legendre knots.
    -> Sampling_points<decltype(curve(start))>
{
    const gsl::Gauss_Legendre g{points};
    Sampling_points<decltype(curve(start))> result(points);
    for (std::size_t i{0}; i<points; ++i) {
        const auto p{g.point(start,end,i)};
        result[i] =
            std::make_tuple(curve(p.first),p.second,derivative(p.first));
    }
    return result;
}

template<typename F1, typename F2>
auto knots_along_piecewise_curve(std::vector<double> boundaries,
        std::vector<std::size_t> points, const F1& curve, const F2& derivative)
    /// @brief Compute `curve` and `derivative` at Gauss-Legendre knots for a
    /// piecewise defined curve.
    ///
    /// @param boundaries cf. `Curve::boundaries()`
    /// @param points The number of knots along the (different segements
    /// of the) curve.
    /// @param curve the curve
    /// @param derivative the derivative of the curve
    -> Sampling_points<decltype(curve(decltype(boundaries)::value_type{}))>
{
    if (boundaries.size() != points.size()+1)
        throw std::invalid_argument{"Each segment requires a number of knots."};
    Sampling_points<decltype(curve(decltype(boundaries)::value_type{}))> result;
    for (std::size_t i{0}; i<points.size(); ++i) {
        auto segment{knots_along_curve(boundaries[i],boundaries[i+1],
                points[i],curve,derivative)};
        result.insert(result.cend(),
                std::make_move_iterator(segment.begin()),
                std::make_move_iterator(segment.end()));
    }
    return result;
}

/// A point in the (x,z)-plane.
struct Point {
    Complex x;
    double x_weight;
    Complex x_derivative;
    double z;
    double z_weight;
};

/// A curve in the complex plane.
struct Curve {
    /// @brief The pair represents the lowest and highest value of a variable
    /// parametrising a curve.
    using Segment = std::optional<std::pair<double,double>>;

    virtual Complex curve_func(double x) const=0;
        ///< Evaluate the curve at `x`.
    virtual Complex derivative_func(double x) const=0;
        ///< Evaluate the derivative of the curve at `x`.
    virtual Segment hits(const Complex& s) const=0;
        ///< @brief Determine whether `s` hits the curve.
        ///<
        ///< If `s` lies on the curve, return the parameter values marking
        ///< the beginning and the end of the segment that is hit.
    virtual std::vector<double> boundaries() const=0;
        ///< @brief Return the boundaries of the parameter values.
        ///<
        ///< In general, the curve is defined piecewisely. `boundaries`
        ///< returns the values of the curve parameters corresponding to
        ///< the start of the curve, the points at which the pieces are glued
        ///< together and the end of the curve. E.g. for a curve connecting
        ///< points A and B in the complex plane, one has
        ///<
        ///<    curve_func(boundaries()[0]) == A
        ///<
        ///<    curve_func(boundaries()[1]) == B
        ///<
        ///< while for a curve connecting A with B and B with C one has
        ///<
        ///<    curve_func(boundaries()[0]) == A
        ///<
        ///<    curve_func(boundaries()[1]) == B
        ///<
        ///<    curve_func(boundaries()[2]) == C
};

template<typename T>
/// A grid in the (x,z)-plane.

/// The z-values are independent of the x-values, that is, for each x-value,
/// the corresponding z-values are the same. The x-values can be specified by
/// an arbitrary curve in the complex plane, while the z-values are straight
/// lines from -1 to 1.
/// This class acts like a decorator for a `Curve` that describes a curve in
/// the x-plane. While the `Curve` is a continuous parametrisation, a `Grid`
/// allows on top of this for discrete sampling with Gauss-Legendre weights.
class Grid : public std::enable_if_t<std::is_base_of<Curve,T>::value,T> {
public:
    Grid(const T& t, std::vector<std::size_t> x_sizes, std::size_t z_size);
        ///< @param t The continuous curve in the x-plane.
        ///< @param x_sizes The number of knots along the (different segements
        ///< of the) curve in the x-plane.
        ///< @param z_size The number of knots along the line in the z-plane.

    Point operator()(std::size_t x_index, std::size_t z_index) const;
        ///< Return the point of the grid at the corresponding position.
    std::vector<double> x_parameter_values() const;
        ///< @brief Return the parameter values at which the curve in the
        ///< x-plane is evaluated at according to the Gauss-Legendre method.
    Complex x(std::size_t x_index) const;
        ///< Return the x-value correspoding to `x_index`.
    Complex derivative(std::size_t x_index) const;
        ///< Return the derivative correspoding to `x_index`.
    double z(std::size_t z_index) const;
        ///< Return the z-value correspoding to `z_index`.
    std::size_t x_size() const noexcept;
        ///< Return the number of knots along the curve in the x-plane.
    std::size_t z_size() const noexcept;
        ///< Return the number of knots along the line in the z-plane.
    double x_parameter_lower() const noexcept;
        ///< @brief Return the parameter corresponding to the beginning of the
        ///< curve in the x-plane.
    double x_parameter_upper() const noexcept;
        ///< @brief Return the parameter corresponding to the end of the
        ///< curve in the x-plane.
private:
    constexpr static double z_lower{-1.0};
    constexpr static double z_upper{1.0};

    double _x_lower;
    double _x_upper;
    std::vector<size_t> x_sizes;
    Sampling_points<Complex> x_knots;
    Knots z_knots;
};

template<typename T>
Grid<T>::Grid(const T& t, std::vector<std::size_t> x_sizes, std::size_t z_size)
    : T{t},
    _x_lower{t.boundaries().front()},
    _x_upper{t.boundaries().back()},
    x_sizes{x_sizes},
    x_knots{
        knots_along_piecewise_curve(
                t.boundaries(),
                x_sizes,
                [t](const auto x) {return t.curve_func(x);},
                [t](const auto x) {return t.derivative_func(x);})},
    z_knots{generate_knots(z_lower,z_upper,z_size)}
{
}

template<typename T>
Point Grid<T>::operator()(std::size_t i, std::size_t j) const
{
    const auto x{x_knots[i]};
    const auto z{z_knots[j]};
    return {std::get<0>(x),std::get<1>(x),std::get<2>(x),z.first,z.second};
}

template<typename T>
std::vector<double> Grid<T>::x_parameter_values() const
{
    const auto identity{[](double x){return x;}};
    auto knots{
        knots_along_piecewise_curve(T::boundaries(),x_sizes,identity,identity)};
    std::vector<double> result(x_size());
    std::transform(knots.cbegin(),knots.cend(),result.begin(),
            [](const auto& t){return std::get<0>(t);});
    return result;
}

template<typename T>
Complex Grid<T>::x(std::size_t x_index) const
{
    return std::get<0>(x_knots[x_index]);
}

template<typename T>
Complex Grid<T>::derivative(std::size_t x_index) const
{
    return std::get<2>(x_knots[x_index]);
}

template<typename T>
double Grid<T>::z(std::size_t z_index) const
{
    return z_knots[z_index].first;
}

template<typename T>
std::size_t Grid<T>::x_size() const noexcept
{
    return x_knots.size();
}

template<typename T>
std::size_t Grid<T>::z_size() const noexcept
{
    return z_knots.size();
}

template<typename T>
double Grid<T>::x_parameter_lower() const noexcept
{
    return _x_lower;
}

template<typename T>
double Grid<T>::x_parameter_upper() const noexcept
{
    return _x_upper;
}

template<typename T>
Grid<T> make_grid(const T& t, std::vector<std::size_t> x_sizes,
        std::size_t z_size)
    /// Construct and return a Grid.
{
    static_assert(std::is_base_of<Curve,T>::value,
            "T needs to inherit from Curve");
    return Grid<T>{t,x_sizes,z_size};
}
} // grid

#endif // KERNEL_GRID_H
