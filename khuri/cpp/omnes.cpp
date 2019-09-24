#include "omnes.h"

namespace omnes {
double derivative_0(const gsl::Function& phase, double threshold, double cut,
        double constant, gsl::Integration& integrate)
    // Return the derivative of the Omnes function at s=0. The parameters are
    // the same as the ones with the same name in the constructor of
    // `class Omnes`.
{
    double first{integrate([&phase](double x){return phase(x)/(x*x);},
            threshold,cut).first};
    double second{constant/cut};
    return (first + second)/constants::pi();
}

Omnes::Omnes(const gsl::Function& phase, double threshold,
        double minimal_distance, gsl::Settings config)
: phase_below{phase},
    constant{0.0}, // value of the `constant` is irrelevant if `cut` is infinity
    threshold{threshold}, cut{std::numeric_limits<double>::infinity()},
    minimal_distance{minimal_distance},
    integrate{config},
    derivative{derivative_0(phase,threshold,cut,constant,integrate)}
{
}

Omnes::Omnes(const gsl::Function& phase, double constant, double threshold,
        double cut, double minimal_distance, gsl::Settings config)
: phase_below{phase}, constant{constant}, threshold{threshold}, cut{cut},
    minimal_distance{minimal_distance},
    integrate{config},
    derivative{derivative_0(phase,threshold,cut,constant,integrate)}
{
}

std::complex<double> Omnes::operator()(const std::complex<double>& s) const
{
    // Apply the Schwartz reflection principle.
    if (s.imag()<0)
        return std::conj(upper(std::conj(s)));
    else
        return upper(s);
}

std::complex<double> Omnes::upper(const std::complex<double>& s) const
{
    if (hits_threshold(s))
        return threshold_presciption(s.real());
    if (hits_cut(s)) 
        return cut_prescription(s.real());
    else
        return ordinary_prescription(s);
}

bool Omnes::hits_threshold(const std::complex<double>& s) const
{
    double distance{std::abs(s-threshold)};
    return distance<=minimal_distance;
}

bool Omnes::hits_cut(const std::complex<double>& s) const
{
    return s.real()>=threshold && std::abs(s.imag())<=minimal_distance;
}

std::complex<double> Omnes::threshold_presciption(double s) const
{
    // average
    return (cut_prescription(threshold+minimal_distance)
            +ordinary_prescription(threshold-minimal_distance)) / 2.0;
}

std::complex<double> Omnes::ordinary_prescription(
        const std::complex<double>& s) const
{
    std::complex<double> above_cut{std::log(1.0-s/cut)};
    auto integral{std::get<0>(cauchy::c_integrate(
                [&s,this](double z){return phase_below(z)/(z*(z-s));},
                threshold,cut,integrate))};
    return std::exp((s*integral-constant*above_cut)/constants::pi());
}

std::complex<double> Omnes::cut_prescription(double s) const
{
    return abs_cut(s)*std::exp(std::complex<double>{0,phase(s)});
}

double Omnes::phase(double s) const
{
    if (s<cut)
        return phase_below(s);
    else
        return constant;
}

double abs_helper(double s, double value)
    // Simplify the calculation of the absolute value of the Omnes
    // function along the cut.
{
    double temp{1.0 - s/value};
    return std::log(std::abs(1.0/temp));
}

double Omnes::abs_cut(double s) const
{
    double phase_at_s{phase(s)};
    auto integral{integrate(
                [&s,&phase_at_s,this](double z)
                    {return (phase_below(z)-phase_at_s)/(z*(z-s));},
                threshold,cut).first};
    double a{s<cut ? constant-phase_at_s : 0.0};
    return std::exp((s*integral + a*abs_helper(s,cut)
            + phase_at_s*abs_helper(s,threshold))/constants::pi());
}
} // omnes
