#include "piecewise.h"

namespace piecewise {
Piecewise::Piecewise(const std::vector<Complex>& knots,
        const std::vector<Para>& parametrisations)
    : parametrisations{parametrisations}
{
    const auto s{parametrisations.size()};
    if (s+1 != knots.size())
        throw std::invalid_argument{
            "Each curve segment needs one parametrisation."};

    std::vector<Complex> differences(knots.size());
    std::adjacent_difference(knots.cbegin(),knots.cend(),differences.begin());

    pieces.resize(s);
    adjacent.resize(s);

    for (std::size_t i{0}; i<s; ++i) {
        pieces[i] = std::make_pair(differences[i+1],knots[i]);
        adjacent[i] = std::make_pair(knots[i],knots[i+1]);
    }
}

std::size_t Piecewise::piece_index(double x) const
{
    if (x<lower() || x>upper())
        throw std::out_of_range{
            "Tried to evaluate piecewise curve outside domain of definition."};
    const auto index{static_cast<std::size_t>(x)};
    // class invariant assures that upper()>=1.
    return index==upper() ? index-1 : index;
}

Complex Piecewise::curve_func(double x) const
{
    const auto k{piece_index(x)};
    const auto& p{pieces.at(k)};
    switch (parametrisations[k]) {
        case linear:
            return p.first*(x-k) + p.second;
        case quadratic:
            return p.first*square(x-k) + p.second;
        default:
            throw Unknown_para{};
    }
}

Complex Piecewise::derivative_func(double x) const
{
    const auto k{piece_index(x)};
    const auto& p{pieces.at(k)};
    switch (parametrisations[k]) {
        case linear:
            return p.first;
        case quadratic:
            return 2.0*p.first*(x-k);
        default:
            throw Unknown_para{};
    }
}

bool in_between(const Complex& x, const Complex& a, const Complex& b)
    // Determine whether `x` is on the connection line of `a` and `b`
    // in between `a` and `b`.
{
    constexpr double minimal_distance{1e-10};
    const double difference{std::abs(x-a) + std::abs(x-b) - std::abs(a-b)};
    return std::abs(difference) < minimal_distance;
}

Piecewise::Segment Piecewise::hits(const Complex& s) const
{
    auto position{std::find_if(adjacent.cbegin(),adjacent.cend(),
            [s](const auto& p){return in_between(s,p.first,p.second);})};
    if (position==adjacent.cend())
        return std::nullopt;
    double lower = std::distance(adjacent.cbegin(),position);
    double upper{lower+1.0};
    return std::make_pair(lower,upper);
}

std::vector<double> Piecewise::boundaries() const
{
    std::vector<double> result(pieces.size()+1);
    std::iota(result.begin(),result.end(),lower());
    return result;
}

std::vector<Complex> vector_decay_points(double pion_mass, double virtuality,
        double cut)
{
    const double m2{square(pion_mass)};
    const double a{virtuality-2.5*m2};
    const double b{-7.0*m2};

    const double x1{4.0*m2};
    const Complex x2{5.0*m2,b};
    const Complex x3{a,b};
    const double x4{a};
    const double x5{mandelstam::s_greater(pion_mass,virtuality)};

    return {x1,x2,x3,x4,x5,cut};
}

Vector_decay::Vector_decay(double pion_mass, double virtuality, double cut)
    : Piecewise{vector_decay_points(pion_mass,virtuality,cut),
        Piecewise::all_linear(5)}
{
}

std::vector<Complex> adaptive_points(double pion_mass, double virtuality,
        double cut)
{
    const auto m2{square(pion_mass)};
    const mandelstam::Critical critical{pion_mass,virtuality};
    const auto lower{-critical.imaginary_radius()};
    const auto right{critical.right()+m2};

    const auto x1{4.0*m2};
    const Complex x2{x1,lower};
    const Complex x3{right,lower};
    const auto x4{right};
    const auto x5{mandelstam::s_greater(pion_mass,virtuality)};

    return {x1,x2,x3,x4,x5,cut};
}

Adaptive::Adaptive(double pion_mass, double virtuality, double cut)
    : Piecewise{adaptive_points(pion_mass,virtuality,cut),
                Piecewise::all_linear(5)}
{
}
} // piecewise
