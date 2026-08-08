/** (c) 2024 Kangrui Xue
 *
 * \file FluidSound.h
 * \brief Public interface for FluidSound library; declares Solver class
 * 
 * Based on code from Timothy Langlois and Ryan Aronson. Thanks to Zhehao Li for reviewing!
 */

#ifndef FLUID_SOUND_H
#define FLUID_SOUND_H

#include <set>

#include "Integrators.h"


namespace FluidSound {

/**
 * \class Solver
 * \brief High-level manager for bubble-based water sound synthesis
 */
template <typename T>
class Solver
{
public:
    /**
     * \brief Constructor: reads Bubble data from file and initializes Oscillators
     * \param[in]  bubFile      path to full bubble tracking file (complete merge/split graph)
     * \param[in]  filteredFile optional path to filtered file; if non-empty, only bubbles in this
     *                          file contribute to output (others still participate in coupling)
     * \param[in]  dt           timestep size
     * \param[in]  scheme       coupling scheme (0 - uncoupled, 1 - coupled)
     * \param[in]  ts           simulation start time (default 0.)
     * \param[in]  timeJitterHalfWidth  if > 0, each oscillator timeline is shifted by U(-w,w) to desynchronize grid artifacts
     * \param[in]  timeJitterSeed       seed for jitter RNG (0 = non-deterministic from std::random_device)
     * \param[in]  denseEvents          if true, every per-sample-line solveTime is inserted into the
     *                                  integrator's event-time set so that K=w0^2 is sampled at every
     *                                  trackedBubInfo row (instead of being a linear ramp between each
     *                                  oscillator's start and end time, which is the default).
     * \param[in]  dampingCoeff         multiplier on the per-sample beta computed by Oscillator::calcBeta.
     *                                  Default 1.0 (original Czerski/Deane radiative+viscous+thermal model);
     *                                  values < 1 lengthen ringdown (longer audible ring), values > 1
     *                                  shorten it.
     * \param[in]  applyListenerAttenuation  if true, multiply each oscillator's contribution to
     *                                  the per-sample sum by 1 / max(|pos(t) - listenerPos|, eps),
     *                                  where pos(t) is the bubble's interpolated trackedBubInfo
     *                                  position at the current time and eps is a small positive
     *                                  clamp (default 1e-6) protecting against listener-on-bubble.
     *                                  Default false (no attenuation; faithful to legacy behavior).
     * \param[in]  listenerX, listenerY, listenerZ  listener position in the same coordinate space
     *                                  as the bubble positions in bubFile. Only used when
     *                                  applyListenerAttenuation is true.
     * \param[in]  listenerEpsilon      lower clamp on distance for the 1/d attenuation, in the
     *                                  same units as bubble positions. Default 1e-6.
     */
    Solver(const std::string& bubFile, const std::string& filteredFile, double dt, int scheme, double ts = 0.,
        double timeJitterHalfWidth = 0., unsigned long long timeJitterSeed = 0ULL,
        double transientPeriods = 0., double transientGain = 1., double forcingCutoff = 0.0006,
        ForcingEnvelope forcingEnvelope = ForcingEnvelope::HARD,
        bool denseEvents = false,
        double dampingCoeff = 1.0,
        unsigned forcingTypeMask = 0x7u,
        const std::string& eventLogPath = std::string(),
        bool applyListenerAttenuation = false,
        double listenerX = 0., double listenerY = 0., double listenerZ = 0.,
        double listenerEpsilon = 1e-6);

    /** \brief Timesteps Oscillator vibrations */
    T step();

    //void loadState(const std::string &stateFile);
    //void saveState(const std::string &stateFile);

    /** \brief Returns vector of ALL Oscillators, sorted by start time */
    std::vector<Oscillator<T>>& oscillators() { return _oscillators; }
    
    /** \brief Returns vector of sorted event times (i.e., when bubbles are added or removed) */
    const std::vector<double>& eventTimes() { return _eventTimes; }

    /** \brief Prints timings from Integrator */
    void printTimings()
    {
        std::cout << "K,C,F time: " << _integrator->coeff_time.count() << std::endl;
        std::cout << "M^-1 time:  " << _integrator->mass_time.count() << std::endl;
        std::cout << "Solve time: " << _integrator->solve_time.count() << std::endl;
    }

    /** \brief Write the --event-log CSV (with measured per-impulse peak |v''|); call after synthesis. */
    void writeEventLog() const;

    ~Solver() { delete _integrator; }

private:
    std::string _eventLogPath;  //!< --event-log destination ("" = disabled)
    double _dt = 0.;    //!< timestep size
    double _ts = 0.;    //!< simulation start time
    int _step = 0;      //!< current time step

    Integrator<T>* _integrator;

    std::vector<Oscillator<T>*> _coupled_osc;
    std::vector<Oscillator<T>*> _uncoupled_osc;

    std::vector<Oscillator<T>> _oscillators;   //!< vector of ALL Oscillators, sorted by start time
    int _osID = 0;  //!< current _oscillators index

    std::vector<double> _eventTimes;    //!< vector of sorted event times (i.e., when to refactor the mass matrix)
    int _evID = 0;  //!< current _eventTimes index

    std::set<int> _contributingBubIDs;   //!< if non-empty, only oscillators with bubIDs in this set contribute to output

    /** \brief If true, each oscillator's per-step accel is scaled by 1 / max(d, eps) with
     *  d = |interp(time).xyz - listenerPos|. Coupling/integrator dynamics are NOT scaled --
     *  only the contribution to the audio sum returned by step(). This gives a clean
     *  far-field listener attenuation without re-coupling the bubble cloud.
     */
    bool _applyListenerAttenuation = false;
    T _listenerX = T(0);
    T _listenerY = T(0);
    T _listenerZ = T(0);
    T _listenerEpsilon = T(1e-6);

    /**
     * \private Given bubble data, chains Bubbles together to form Oscillators
     * \param[in]  bubMap  map from Bubble IDs to Bubble objects
     */
    void _makeOscillators(const std::map<int, Bubble<T>>& bubMap, double timeJitterHalfWidth, unsigned long long timeJitterSeed,
        double transientPeriods, double transientGain, double forcingCutoff, bool denseEvents, double dampingCoeff,
        unsigned forcingTypeMask, const std::string& eventLogPath);
};

} // namespace FluidSound

#endif // #ifndef FLUID_SOUND_H