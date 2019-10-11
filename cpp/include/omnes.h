#ifndef OMNES_FUNCTION_H
#define OMNES_FUNCTION_H

#include "cauchy.h"
#include "constants.h"
#include "gsl_interface.h"

#include <cmath>
#include <complex>

/// The Omnes function of an arbitrary phase.

/// The way of thinking about the Omnes function is that it provides a
/// different function for different phases.
/// Hence, a class `Omnes` is provided. Each instance of the class provides
/// a different Omnes function, i.e. one for a different specific phase.
namespace omnes {
/// The Omnes function for arbitrary phases and thresholds.
class Omnes {
public:
    Omnes(const gsl::Function& phase, double threshold,
            double minimal_distance, gsl::Settings config=gsl::Settings{});
        ///< @param phase The phase of the Omnes function above its branch cut.
        ///< @param threshold The start of the branch cut of the Omnes function.
        ///< @param minimal_distance Half the width of a band
        ///< around the cut. For arguments of the Omnes function in this band,
        ///< a different prescription is used for the evaluation of the Omnes
        ///< function to take care of the singularity in the integral.
        ///< @param config The settings for the integration routine.

    Omnes(const gsl::Function& phase, double threshold, double constant,
            double cut, double minimal_distance,
            gsl::Settings config=gsl::Settings{});
        ///< @param phase The phase of the Omnes function in
        ///< [`threshold`,`cut`].
        ///< @param threshold The start of the branch cut of the Omnes function.
        ///< @param constant The phase of the Omnes function is set equal to
        ///< `constant` along the real line above `cut`.
        ///< @param cut Cf. `phase` and `constant`
        ///< @param minimal_distance Half the width of a band
        ///< around the cut. For arguments of the Omnes function in this band,
        ///< a different prescription is used for the evaluation of the Omnes
        ///< function to take care of the singularity in the integral.
        ///< @param config The settings for the integration routine.

    std::complex<double> operator()(std::complex<double> s) const;
        ///< Evaluate the Omnes function at `s`.

    double derivative_at_zero() const noexcept {return derivative;}
        ///< Return the derivative of the Omnes function at the origin.
private:
    const gsl::Function phase_below; // phase below `cut`
    const double constant;
    const double threshold;
    const double cut;
    const double minimal_distance;
    const gsl::Cquad integrate;
    const double derivative;

    std::complex<double> upper(const std::complex<double>& s) const;
        // Evaluate the Omnes function in the upper half of the complex plane
    bool hits_threshold(const std::complex<double>& s) const;
        // Return true if `s` is close to the `threshold`, false otherwise.
    bool hits_cut(const std::complex<double>& s) const;
        // Return true if `s` is in the region around the branch cut, false
        // otherwise.
    std::complex<double> threshold_presciption(double s) const;
        // Calculate the Omnes function if `s` is close to `threshold`.
    std::complex<double>
            ordinary_prescription(const std::complex<double>& s) const;
        // Calculate the Omnes function if `s` is not close to the branch
        // cut.
    std::complex<double> cut_prescription(double s) const;
        // Calculate the Omnes function if `s` is close to the branch cut.
    double phase(double s) const;
        // Calculate the phase of the Omnes function along the branch cut.
    double abs_cut(double s) const;
        // Calculate the absolute value of the Omnes function along the branch
        // cut.
};
} // omnes

#endif // OMNES_FUNCTION_H
