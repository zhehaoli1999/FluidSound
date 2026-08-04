/** (c) 2024 Kangrui Xue
 *
 * \file Oscillator.cpp
 * \brief Implements Oscillator class forcing and damping functions
 * 
 * References:
 *   [Langlois et al. 2016] Toward Animating Water with Complex Acoustic Bubbles
 */

#include <algorithm>
#include <cmath>
#include <random>

#include "Oscillator.h"


namespace FluidSound {

static const double RHO_WATER = 998.;	// density of water
static const double GAMMA = 1.4;		// gas heat capacity ratio
static const double SIGMA = 0.0726;		// surface tension
static const double MU = 8.9e-4;		// dynamic viscosity of water
static const double GTH = 1.6e6;		// thermal damping constant
static const double G = 1.0;			// TODO: ?
static const double CF = 1497;			// speed of sound in water
static const double ATM = 101325;		// atmospheric pressure


static std::default_random_engine s_forcingRnd;
static std::uniform_real_distribution<double> s_eta(0.4, 1.5);
static std::uniform_real_distribution<double> s_frac(0.4, 0.8);


/** Eq. 14 from [Langlois et al. 2016] */
template <typename T>
std::pair<T, T> Oscillator<T>::CzerskiJetForcing(T radius, T maxCutoff)
{
    // T eta = 0.95;  // TODO
    T eta = 0.95;  // TODO

    T cutoff = std::min(maxCutoff, T(0.5) / (T(3) / radius));   // 1/2 minnaert period

    T pressure_in0 = (ATM + 2. * SIGMA / radius);
    T weight = -9. * GAMMA * SIGMA * eta * pressure_in0 * std::sqrt(1. + eta * eta) / (4. * RHO_WATER * radius * radius * radius);

    T mass = (RHO_WATER / (4. * M_PI * radius));

    return std::pair<T, T>(cutoff, weight / mass);
}

/** Eq. 15 from [Langlois et al. 2016] */
template <typename T>
std::pair<T, T> Oscillator<T>::MergeForcing(T radius, T r1, T r2, T maxCutoff)
{
    T frac = s_frac(s_forcingRnd);
    T factor = std::pow(2. * SIGMA * r1 * r2 / (RHO_WATER * (r1 + r2)), 0.25);

    T cutoff = std::min(maxCutoff, T(0.5) / (T(3) / radius));   // 1/2 minnaert period
    T tmp = std::pow(frac * std::min(r1, r2) / 2. / factor, 2);     // TODO: cleanup
    cutoff = std::min(cutoff, tmp);

    T pressure_in0 = (ATM + 2. * SIGMA / radius);
    T weight = 6. * SIGMA * GAMMA * pressure_in0 / (RHO_WATER * radius * radius * radius);

    T mass = (RHO_WATER / (4. * M_PI * radius));

    return std::pair<T, T>(cutoff, weight / mass);
}

/** See Oscillator.h. Unit-weight response integration: with F(t)/m =
 *  envelope(t/tau) * t^2 (weight 1), the response v_u is linear in weight, so
 *  v_peak(w) = w * v_peak_u and the injected work W(w) = w^2 * m * I with
 *  I = int envelope * t^2 * v_u' dt. One scalar RK4 per event at load time. */
template <typename T>
std::pair<T, T> Oscillator<T>::CalibratedForcing(T radius, T w0, T beta, T weightSign,
    T deltaArea, const ForcingParams& fp, bool smoothstepEnvelope)
{
    double r = static_cast<double>(radius);
    double w = static_cast<double>(w0);
    double b = static_cast<double>(beta);
    if (!(r > 0.) || !(w > 0.)) { return std::pair<T, T>(T(0), T(0)); }

    double period = 2. * M_PI / w;
    double tau = std::min(std::max(fp.zeta * period, fp.tauMin), fp.tauMax);

    double eps0 = fp.epsRef * std::pow(r / fp.epsRRef, -fp.epsExp);
    eps0 = std::min(std::max(eps0, fp.epsMin), fp.epsMax);
    double vTarget = 4. * M_PI * r * r * r * eps0;   // dV = 3*eps0*V0

    // Unit-weight response over the forcing window + two ringdown periods.
    double Tend = tau + 2. * period;
    int steps = std::max(4000, static_cast<int>(Tend / period * 64.));
    steps = std::min(steps, 400000);
    double h = Tend / steps;

    auto Fm = [&](double t) -> double {
        if (t < 0. || t >= tau) { return 0.; }
        double envelope = 1.;
        if (smoothstepEnvelope && tau > 0.)
        {
            double x = t / tau;
            envelope = 1. - x * x * (3. - 2. * x);
        }
        return envelope * t * t;
    };
    auto acc = [&](double t, double v, double vd) -> double {
        return Fm(t) - 2. * b * vd - w * w * v;
    };

    double v = 0., vd = 0., vPeak = 0., I = 0., prevPow = 0.;
    for (int i = 0; i < steps; i++)
    {
        double t = i * h;
        double k1v = vd,                 k1a = acc(t, v, vd);
        double k2v = vd + 0.5 * h * k1a, k2a = acc(t + 0.5 * h, v + 0.5 * h * k1v, vd + 0.5 * h * k1a);
        double k3v = vd + 0.5 * h * k2a, k3a = acc(t + 0.5 * h, v + 0.5 * h * k2v, vd + 0.5 * h * k2a);
        double k4v = vd + h * k3a,       k4a = acc(t + h, v + h * k3v, vd + h * k3a);
        v  += h / 6. * (k1v + 2. * k2v + 2. * k3v + k4v);
        vd += h / 6. * (k1a + 2. * k2a + 2. * k3a + k4a);

        vPeak = std::max(vPeak, std::abs(v));
        double pow_i = Fm(t + h) * vd;               // unit-weight F/m times v'
        I += 0.5 * (prevPow + pow_i) * h;
        prevPow = pow_i;
    }
    if (!(vPeak > 0.)) { return std::pair<T, T>(T(0), T(0)); }

    double weight = vTarget / vPeak;

    // Surface-tension energy bound (Doug's falsification bound, enforced):
    // an event that changes the interface area by dA released at most
    // sigma*|dA|; the breathing mode may receive only a fraction of it.
    if (fp.energyCapEta > 0.)
    {
        double dA = (deltaArea > 0.) ? static_cast<double>(deltaArea) : 4. * M_PI * r * r;
        double mass = RHO_WATER / (4. * M_PI * r);
        double Epred = weight * weight * mass * I;
        double Ecap = fp.energyCapEta * SIGMA * dA;
        if (Epred > Ecap && Epred > 0.)
        {
            weight *= std::sqrt(Ecap / Epred);
        }
    }

    return std::pair<T, T>(static_cast<T>(tau), static_cast<T>(weightSign) * static_cast<T>(weight));
}

/**  */
template <typename T>
T Oscillator<T>::calcBeta(T radius, T w0, T coeff)
{
    T dr = w0 * radius / CF;
    T dvis = 4 * MU / (RHO_WATER * w0 * radius * radius);
    T phi = 16. * GTH * G / (9 * (GAMMA - 1) * (GAMMA - 1) * w0 / 2. / M_PI);
    T dth = 2 * (std::sqrt(phi - 3) - (3 * GAMMA - 1) / (3 * (GAMMA - 1))) / (phi - 4);

    T dtotal = dr + dvis + dth;

    return coeff * w0 * dtotal / std::sqrt(dtotal * dtotal + 4);
}

template class Oscillator<float>;
template struct Oscillator<double>;

}