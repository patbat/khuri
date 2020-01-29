#include "mandelstam.h"

namespace mandelstam {
double unit(double pion_mass, double virtuality)
    // helper function for parametrization of complex egg
{
    static double sqrt2{std::sqrt(2.0)};
    return sqrt2
        * std::sqrt(s_greater(pion_mass,virtuality)-4.0*pion_mass*pion_mass);
}

double change_1(double x, double pion_mass) noexcept
    // change of variables in parametrization of complex egg
{
    return 4.0*pion_mass*pion_mass + x*x/4.0;
}

double change_2(double x, double pion_mass, double virtuality)
    // change of variables in parametrization of complex egg
{
    double temp{2*unit(pion_mass,virtuality)-x};
    return s_greater(pion_mass,virtuality) - temp*temp/4.0;
}

void inside_region(double x, double unit)
{
    if (x < 0.0 || 2.0*unit < x)
        throw std::invalid_argument{"Egg is not defined in this region."};
}

Egg::Egg(double pion_mass, double virtuality)
    : pion_mass{pion_mass}, virtuality{virtuality},
    s_greater_{s_greater(pion_mass,virtuality)},
    s_smaller_{s_smaller(pion_mass,virtuality)},
    unit_{unit(pion_mass,virtuality)}
{
}

Complex Egg::lower_segment(double x) const
{
    inside_region(x,unit_);
    double y{
        x<=unit_ ? change_1(x,pion_mass) : change_2(x,pion_mass,virtuality)};
    return t_photon_pion_min(y,pion_mass,virtuality);
}

Complex Egg::upper_segment(double x) const
{
    return std::conj(lower_segment(4.0*unit_ - x));
}

Complex Egg::operator()(double x) const
{
    if (x <= change())
        return lower_segment(x);
    return upper_segment(x);
}

Complex Egg::first_half(double x) const
{
    double y{change_1(x,pion_mass)};
    double sig{sigma(y,pion_mass).real()};
    double sq{std::sqrt((y-s_smaller_)*(s_greater_-y))};
    double real{-x/4.0};
    double imag{pion_mass*pion_mass/(y*y)*sq*std::sqrt(y)
        + x/8.0*sig*(s_greater_+s_smaller_-2.0*y)/sq};
    return {real,-imag};
}

Complex Egg::second_half(double x) const
{
    double y{change_2(x,pion_mass,virtuality)};
    double sig{sigma(y,pion_mass).real()};
    double sq{std::sqrt(y-s_smaller_)};
    double shift{x/2.0-unit_};
    double real{shift/2.0};
    double m2{pion_mass*pion_mass}; 
    double imag{-shift*m2/(y*y)*sq*std::sqrt(y*(s_greater_-y)/(y-4.0*m2))
        + sig/4.0*(s_greater_+s_smaller_-2.0*y)/sq};
    return {real,-imag};
}

Complex Egg::lower_derivative(double x) const
{
    inside_region(x,unit_);
    if (x <= unit_)
        return first_half(x);
    return second_half(x);
}

Complex Egg::upper_derivative(double x) const
{
    return -std::conj(lower_derivative(4.0*unit_ - x));
}

Complex Egg::derivative(double x) const
{
    if (x <= change())
        return lower_derivative(x);
    return upper_derivative(x);
}

double Egg::lower(double s) const
{
    double threshold{4*pion_mass*pion_mass};
    if (s < threshold || s_greater_ < s) 
        throw std::invalid_argument{"Egg is not defined in this region."};
    double boundary{(s_greater_+threshold) / 2.0};
    if (s<boundary)
        return 2.0*std::sqrt(s-threshold);
    return 2.0 * (unit_ - std::sqrt(s_greater_-s));
}

double Egg::upper(double s) const
{
    return 4.0*unit_ - lower(s);
}
} // mandelstam
