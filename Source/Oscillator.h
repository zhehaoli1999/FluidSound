/** (c) 2024 Kangrui Xue
 *
 * \file Oscillator.h
 * \brief Declares Oscillator struct, including computation of forcing and damping terms
 */

#ifndef _FS_OSCILLATOR_H
#define _FS_OSCILLATOR_H

#include <Eigen/Dense>

#include "BubbleUtils.h"


namespace FluidSound {

/**
 * \class Oscillator
 * \brief Represents a single oscillator, \f$ \ddot{v}(t) + 2\beta \dot{v}(t) + \omega_0^2 v(t) = p(t) / m \f$
 * 
 * Whereas the Bubble struct corresponds to physical bubbles, the Oscillator struct is more of a
 *   mathematical abstraction, meant to interface efficiently with the Integrator class
 */
template <typename T>
struct Oscillator
{
    /** \brief vector of IDs of Bubbles belonging to this Oscillator, sorted by increasing start time */
    std::vector<int> bubIDs;

    double startTime = -1.;
    double endTime = -1.;

    /** \brief current volume displacement and volume velocity state vector: [v v'] */
    Eigen::Vector2<T> state = { 0., 0. };
    T accel = 0.;   /**< \brief current volume acceleration : v'' */

    std::vector<double> solveTimes;  //!< times [ t(0) ... t(N) ] corresponding to solveData
    Eigen::Array<T, 6, Eigen::Dynamic> solveData;
    /**< \brief For indices (0, ..., N), solveData is given by:
     * [[ r(0)  ... r(N)  ]
     *  [ ω0(0) ... ω0(N) ]
     *  [ x(0)  ... x(N)  ]
     *  [ y(0)  ... y(N)  ]
     *  [ z(0)  ... z(N)  ]
     *  [ 2β(0) ... 2β(N) ]]
     */

    /** \brief Returns array of linearly interpolated solve data at specified time */
    Eigen::Array<T, 6, 1> interp(double time)
    {
        if (time >= solveTimes.back()) { return solveData.col(solveTimes.size() - 1); }
        else if (time <= solveTimes[0]) { return solveData.col(0); }

        while (time < solveTimes[_idx]) { _idx--; }
        while (_idx < solveTimes.size() - 1 && time > solveTimes[_idx + 1]) { _idx++; }

        double alpha = (time - solveTimes[_idx]) / (solveTimes[_idx + 1] - solveTimes[_idx]);
        return (1. - alpha) * solveData.col(_idx) + alpha * solveData.col(_idx + 1);
    }

    /** \brief Returns true if this Oscillator has decayed below the audible floor.
     *
     * The audio sample is \f$ \sum_i \mathrm{osc}_i.\mathrm{accel} \f$ ;
     * removing an oscillator from the active set produces a sample-rate step
     * of size \f$ |\mathrm{accel}| \f$ . For a free (uncoupled, unforced)
     * oscillator,
     * \f[
     *     \mathrm{accel} \;=\; -2\beta\,\dot v \;-\; \omega_0^2\, v
     *                    \;\approx\; -\omega_0^2\, v
     *     \qquad (\beta \ll \omega_0)
     * \f]
     * so the worst-case \f$ |\mathrm{accel}| \f$ over the next period is
     * \f$ \omega_0^2 \cdot \mathrm{envelope} \f$ , where
     * \f$ \mathrm{envelope} = \sqrt{v^2 + (\dot v / \omega_0)^2} \f$
     * is the slowly-varying complex amplitude (and is unbiased by the
     * \f$ v \f$ zero-crossing, unlike \f$ |\mathrm{accel}| \f$ itself).
     *
     * The old threshold \f$ \|\mathrm{state}\| < 10^{-10} \f$ is
     * \f$ K \f$ -blind: for an \f$ \omega_0 = 2\pi \cdot 5\,\text{kHz} \f$
     * bubble it still admits \f$ |\mathrm{accel}| \lesssim 10^9 \cdot 10^{-10}
     * = 0.1 \f$ , producing an audible death-click the instant the oscillator
     * leaves the sum in FluidSound::Solver::step. We now require the
     * worst-case \f$ |\mathrm{accel}| \f$ itself to be below an inaudible floor.
     */
    bool is_dead() const
    {
        const T accel_floor = T(1e-6);
        if (solveData.cols() == 0) { return state.norm() < T(1e-10); }
        T w0 = solveData.row(1).maxCoeff();
        if (!(w0 > T(0))) { return state.norm() < T(1e-10); }
        T envelope = std::sqrt(state(0) * state(0) + (state(1) / w0) * (state(1) / w0));
        return (w0 * w0) * envelope < accel_floor;
    }
    
    bool operator < (const Oscillator<T>& osc) const { return startTime < osc.startTime; }

    /** \brief Fade window (in seconds) for the coupling smoothstep on both sides
     *  of the oscillator's lifetime. Multiplied into off-diagonal mass-matrix
     *  entries by Coupled_Direct::_constructMass so that adding or removing this
     *  oscillator from the coupled set is C^1-smooth instead of a step change
     *  (which would otherwise produce a Schur-complement click in every other
     *  coupled bubble's v''). 1 ms = 48 samples at 48 kHz, short enough that
     *  the temporary "weakly coupled" transient is shorter than a Minnaert
     *  period at all audible frequencies. Set to 0 to restore the original
     *  hard-add / hard-remove behavior.
     */
    static constexpr double COUPLING_FADE = 0.001;

    /** \brief Coupling weight in [0, 1] for this oscillator at the given time.
     *
     *  Hermite smoothstep ramp 0 -> 1 over [startTime, startTime + COUPLING_FADE]
     *  and 1 -> 0 over [endTime - COUPLING_FADE, endTime]. The product
     *  \f$ \alpha_i \alpha_j \f$ multiplies every off-diagonal entry of \f$ M \f$
     *  in Coupled_Direct::_constructMass (diagonal stays at 1). At \f$ \alpha = 0 \f$
     *  the oscillator's row/col of \f$ M \f$ is exactly \f$ e_i \f$, so adding or
     *  removing it produces zero Schur correction on the surviving block of
     *  \f$ M^{-1} \f$ -> no click in any other bubble's \f$ \ddot v \f$.
     *
     *  For oscillators whose lifetime is shorter than COUPLING_FADE the fade-in
     *  and fade-out windows overlap and \f$ \alpha \f$ never reaches 1; physically
     *  defensible since sub-fade-window transient bubbles barely couple to
     *  anything anyway.
     */
    T coupling_alpha(double time) const
    {
        if (COUPLING_FADE <= 0.0) { return T(1); }
        double fin;
        if (time <= startTime) { fin = 0.0; }
        else if (time >= startTime + COUPLING_FADE) { fin = 1.0; }
        else
        {
            double x = (time - startTime) / COUPLING_FADE;
            fin = x * x * (3.0 - 2.0 * x);
        }
        double fout;
        if (time <= endTime - COUPLING_FADE) { fout = 1.0; }
        else if (time >= endTime) { fout = 0.0; }
        else
        {
            double x = (endTime - time) / COUPLING_FADE;
            fout = x * x * (3.0 - 2.0 * x);
        }
        return T(fin * fout);
    }

    Eigen::Array<T, 3, Eigen::Dynamic> forceData;
    /**< \brief For indices (0, ..., F), forceData is given by:
     * [[ forceTime(0) ... forceTime(F) ]
     *  [ cutoff(0)    ... cutoff(F)    ]
     *  [ weight(0)    ... weight(F)    ]]
     *
     * All forcing functions have the form F(t) = (t < cutoff) * weight * t * t
     *  (where t is relative to the force start time : t = time - forceTime)
     */

    std::vector<int> forceTypes;   //!< EventType of each forcing column (index-aligned with forceData and bubIDs)
    std::vector<T> forceRadii;     //!< equilibrium radius of the bubble behind each forcing column
    std::vector<T> forcePeakAccel; //!< measured max |v''| attributed to each forcing column (filled during synthesis)
    int forceCursor = 0;           //!< index of the forcing column currently active during synthesis

    /** 
     * \brief Neck collapse forcing model from Czerksi/Deane [2008; 2010]
     * \param[in]  radius  bubble equilibrium radius
     * \return  (cutoff, weight) pair
     */
    static std::pair<T, T> CzerskiJetForcing(T radius, T maxCutoff = T(0.0006));
    
    /** 
     * \brief Neck expansion forcing model from [Czerski 2011]
     * \param[in]  radius  bubble equilibrium radius
     * \param[in]  r1, r2  radii of parent bubbles
     * \return  (cutoff, weight) pair
     */
    static std::pair<T, T> MergeForcing(T radius, T r1, T r2, T maxCutoff = T(0.0006));
    
    /** \brief Damping via radiative, viscous, and thermal effects.
     *  \param[in]  coeff   multiplier on the returned beta (default 1, the
     *                      original Czerski/Deane value). Use coeff < 1 to
     *                      lengthen ringdown, coeff > 1 to shorten it.
     */
    static T calcBeta(T radius, T w0, T coeff = T(1));
    
private:
    int _idx = 0;
};

} // namespace FluidSound

#endif // _FS_OSCILLATOR_H