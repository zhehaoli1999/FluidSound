/** (c) 2024 Kangrui Xue
 *
 * \file FluidSound.cpp
 * \brief Implements Solver class
 */

#include "FluidSound.h"
#include "BubbleUtils.h"

#include <iomanip>
#include <limits>
#include <random>

namespace FluidSound {

/** */
template <typename T>
Solver<T>::Solver(const std::string &bubFile, const std::string &filteredFile, double dt, int scheme, double ts,
    double timeJitterHalfWidth, unsigned long long timeJitterSeed, double transientPeriods, double transientGain,
    double forcingCutoff, ForcingEnvelope forcingEnvelope, bool denseEvents, double dampingCoeff)
    : _dt(dt), _ts(ts)
{
    std::map<int, Bubble<T>> allBubbles;

    std::cout << "Reading bubble data from \"" << bubFile << "\"" << std::endl;
    BubbleUtils<T>::loadBubbleFile(allBubbles, bubFile);

    if (!filteredFile.empty())
    {
        std::cout << "Reading contributing bubble IDs from \"" << filteredFile << "\"" << std::endl;
        _contributingBubIDs = parseBubbleIDsFromFile(filteredFile);
        size_t totalBubbles = allBubbles.size();
        size_t contributing = _contributingBubIDs.size();
        size_t filteredOut = (contributing <= totalBubbles) ? (totalBubbles - contributing) : 0;
        double pctFiltered = (totalBubbles > 0) ? (100.0 * filteredOut / totalBubbles) : 0.0;
        std::cout << "  " << contributing << " bubbles contribute, " << filteredOut
                  << " filtered out (" << std::fixed << std::setprecision(1) << pctFiltered << "%)" << std::endl;
    }
     
    _makeOscillators(allBubbles, timeJitterHalfWidth, timeJitterSeed, transientPeriods, transientGain, forcingCutoff, denseEvents, dampingCoeff);
    std::cout << "Total number of oscillators = " << _oscillators.size() << std::endl;
    std::cout << "Total number of event times = " << _eventTimes.size()
              << (denseEvents ? "  (dense: per-sample-line K=w0^2)" : "  (sparse: linear K=w0^2 ramp per oscillator lifetime)")
              << std::endl;
    if (timeJitterHalfWidth > 0.)
    {
        std::cout << "Per-oscillator time jitter: uniform on [-" << timeJitterHalfWidth << ", " << timeJitterHalfWidth << "] s";
        if (timeJitterSeed != 0ULL)
            std::cout << " (seed=" << timeJitterSeed << ")";
        std::cout << std::endl;
    }

    if (!_contributingBubIDs.empty())
    {
        size_t contributingOsc = 0;
        for (const auto& osc : _oscillators)
        {
            for (int bid : osc.bubIDs)
            {
                if (_contributingBubIDs.count(bid))
                {
                    contributingOsc++;
                    break;
                }
            }
        }
        size_t totalOsc = _oscillators.size();
        size_t filteredOsc = totalOsc - contributingOsc;
        double pctFilteredOsc = (totalOsc > 0) ? (100.0 * filteredOsc / totalOsc) : 0.0;
        std::cout << "  " << contributingOsc << " oscillators contribute, " << filteredOsc
                  << " oscillators filtered out (" << std::fixed << std::setprecision(1) << pctFilteredOsc << "%)" << std::endl;
    }

    if (_eventTimes.empty())
    {
        throw std::runtime_error("No oscillators created (all bubbles filtered out). Check: solve data, frequency < 18 kHz, duration > 3 periods.");
    }

    switch (scheme)
    {
        case 1: _integrator = new Coupled_Direct<T>(dt, forcingEnvelope); break;
        default: _integrator = new Uncoupled<T>(dt, forcingEnvelope); break;
    }
}

/** */
template <typename T>
T Solver<T>::step()
{
    double time = _dt * _step + _ts;

    // Check if any events (e.g., Bubbles added or removed) will occurr during this timestep
    while (_evID < _eventTimes.size() && time >= _eventTimes[_evID])
    {
        if (time < _eventTimes[_evID + 1])
        {
            double time1 = _eventTimes[_evID]; double time2 = _eventTimes[_evID + 1];

            // Check if any Oscillators have ended by time1. We uncouple them from the bubble cloud but
            //    continue timestepping their oscillations (until they die out) to avoid discontinuities.
            int coupled_idx = 0;
            for (Oscillator<T>* osc : _coupled_osc)
            {
                if (time1 >= osc->endTime)
                {
                    _uncoupled_osc.push_back(osc);
                    continue;
                }
                _coupled_osc[coupled_idx] = osc;
                coupled_idx++;
            }
            _coupled_osc.resize(coupled_idx);

            // Of uncoupled Oscillators, check if any have decayed sufficiently by time1.
            //   These Oscillators are fully removed.
            int uncoupled_idx = 0;
            for (Oscillator<T>* osc : _uncoupled_osc)
            {
                if (osc->is_dead())
                {
                    continue;
                }
                _uncoupled_osc[uncoupled_idx] = osc;
                uncoupled_idx++;
            }
            _uncoupled_osc.resize(uncoupled_idx);

            // Check if any Oscillators will start between time1 and time2.
            //    NOTE: _oscillators sorted by increasing startTime.
            while (_osID < _oscillators.size() && time >= _oscillators[_osID].startTime && 
                   _oscillators[_osID].startTime < time2)
            {
                Oscillator<T>* osc = &(_oscillators[_osID]);
                if (time1 < osc->endTime)
                {
                    _coupled_osc.push_back(osc);
                }
                _osID++;
            }

            // Prepare _integrator for timestepping: transfer data + refactor mass matrix 
            _integrator->updateData(_coupled_osc, _uncoupled_osc, time1, time2);
            _integrator->refactor();
        }
        _evID++;
    }
    _step++;

    std::vector<Oscillator<T>*> total_osc(_coupled_osc.begin(), _coupled_osc.end());
    total_osc.insert(total_osc.end(), _uncoupled_osc.begin(), _uncoupled_osc.end());
    size_t N_total = total_osc.size();

    if (N_total == 0) { return 0.; }
    //if (N_total > 1024) { throw std::runtime_error("Too many bubbles for coupling!"); }


    _integrator->step(time);

    // Unpack _integrator->States() to update Oscillator states accordingly
    // If _contributingBubIDs is set, only add accel from oscillators whose bubIDs are in that set
    T total_response = 0.0;
    for (int i = 0; i < N_total; i++)
    {
        total_osc[i]->state(0) = _integrator->States()(i);
        total_osc[i]->state(1) = _integrator->States()(i + N_total);

        total_osc[i]->accel = _integrator->Derivs()(i + N_total);
        if (_contributingBubIDs.empty())
        {
            total_response += total_osc[i]->accel;
        }
        else
        {
            for (int bid : total_osc[i]->bubIDs)
            {
                if (_contributingBubIDs.count(bid))
                {
                    total_response += total_osc[i]->accel;
                    break;
                }
            }
        }
    }
    if (std::abs(total_response) > 100.) { throw std::runtime_error("Instability detected!"); }

    return total_response;
}

/** */
template <typename T>
void Solver<T>::_makeOscillators(const std::map<int, Bubble<T>> &bubMap, double timeJitterHalfWidth, unsigned long long timeJitterSeed,
    double transientPeriods, double transientGain, double forcingCutoff, bool denseEvents, double dampingCoeff)
{
    _oscillators.clear();

    std::set<int> usedBubIDs;
    size_t transientOscillators = 0;
    for (const std::pair<int, Bubble<T>>& bubPair : bubMap)  // Loop over all bubbles
    {
        // Skip if bubble has already been used
        if (usedBubIDs.count(bubPair.first)) continue;
        
        int curBubID = bubPair.first;
        const Bubble<T>* curBub = &bubPair.second;

        // Initialize Oscillator - we will set its data as we go 
        Oscillator<T> osc;
        osc.startTime = curBub->startTime;

        // Temporary buffers for solve data and force data
        std::vector<double> solveTimes;
        std::vector<T> s_radii, s_w0, s_x, s_y, s_z;
        std::vector<T> Cvals;   // precomputed damping coefficient
        std::vector<T> forceTimes, f_cutoffs, f_weights;

        while (true)
        {
            bool lastBub = true;
            usedBubIDs.insert(curBubID);
            
            // Add the solve data for this bubble (if there is any)
            solveTimes.insert(solveTimes.end(), curBub->solveTimes.begin(), curBub->solveTimes.end());
            s_radii.insert(s_radii.end(), curBub->solveTimes.size(), curBub->radius);
            s_w0.insert(s_w0.end(), curBub->w0.begin(), curBub->w0.end());
            s_x.insert(s_x.end(), curBub->x.begin(), curBub->x.end());
            s_y.insert(s_y.end(), curBub->y.begin(), curBub->y.end());
            s_z.insert(s_z.end(), curBub->z.begin(), curBub->z.end());

            osc.bubIDs.push_back(curBubID);
            
            // ----- Handle bubble start event: forcing logic -----
            std::pair<T, T> force(T(0), T(0));

            if (curBub->startType == EventType::SPLIT)
            {
                // Tracker glitches in the LBM exhale dataset frequently produce
                // microscopic-parent -> macroscopic-child SPLIT events (e.g. a 20um
                // ghost bubble "splits" into a 10mm child). The original guard
                // `parent.radius >= child.radius` enforced mass-conservation and
                // therefore zeroed forcing on these glitches, which silenced the
                // child entirely. We now treat every SPLIT child as a freshly
                // entrained bubble and force it via CzerskiJetForcing regardless
                // of the parent's recorded radius.
                //int parentBubID = curBub->prevBubIDs.at(0);
                //if (bubMap.at(parentBubID).radius >= curBub->radius)
                //{
                    force = Oscillator<T>::CzerskiJetForcing(curBub->radius, T(forcingCutoff));
                //}
            }
            else if (curBub->startType == EventType::MERGE)  // TODO: cleanup this code
            {
                if (curBub->prevBubIDs.size() == 2)
                {
                    bool allMerge = true;
                    for (int parentBubID : curBub->prevBubIDs)
                    {
                        allMerge = allMerge && (bubMap.at(parentBubID).endType == EventType::MERGE);
                    }

                    if (allMerge)
                    {
                        int p1 = curBub->prevBubIDs.at(0);
                        int p2 = curBub->prevBubIDs.at(1);

                        T r1 = bubMap.at(p1).radius;
                        T r2 = bubMap.at(p2).radius;

                        if (r1 + r2 > curBub->radius)
                        {
                            T v1 = 4. / 3. * M_PI * r1 * r1 * r1;
                            T v2 = 4. / 3. * M_PI * r2 * r2 * r2;
                            T vn = 4. / 3. * M_PI * curBub->radius * curBub->radius * curBub->radius;

                            T diff = v1 + v2 - vn;
                            if (diff <= std::max(v1, v2))
                            {
                                if (v1 > v2) { v1 -= diff; }
                                else { v2 -= diff; }

                                r1 = std::pow(3. / 4. / M_PI * v1, 1. / 3.);
                                r2 = std::pow(3. / 4. / M_PI * v2, 1. / 3.);

                                force = Oscillator<T>::MergeForcing(curBub->radius, r1, r2, T(forcingCutoff));
                            }
                        }
                        else
                        {
                            force = Oscillator<T>::MergeForcing(curBub->radius, r1, r2, T(forcingCutoff));
                        }
                    }
                }
            }
            else if (curBub->startType == EventType::ENTRAIN)
            {
                force = Oscillator<T>::CzerskiJetForcing(curBub->radius, T(forcingCutoff));
            }
            forceTimes.push_back(curBub->startTime);
            f_cutoffs.push_back(force.first);
            f_weights.push_back(force.second);


            // ----- Handle bubble end event: chaining logic -----
            if (curBub->endType == EventType::MERGE)
            {
                const Bubble<T>& nextBub = bubMap.at(curBub->nextBubIDs.at(0));
                int largestParent = BubbleUtils<T>::largestBubbleID(nextBub.prevBubIDs, bubMap);
                
                if (largestParent == curBubID)
                {
                    lastBub = false;
                    curBubID = curBub->nextBubIDs.at(0);
                    curBub = &bubMap.at(curBubID);
                }
            }
            else if (curBub->endType == EventType::SPLIT)
            {
                // Sort children in order of size
                std::multimap<T, int, std::greater<T>> children;
                for (int childID : curBub->nextBubIDs)
                {
                    children.insert(std::make_pair(bubMap.at(childID).radius, childID));
                }
                // Continue this bubble to the largest child bubble where this is the largest parent
                for (auto &child : children)
                {
                    if (usedBubIDs.count(child.second)) continue;
                    
                    int largestParent = BubbleUtils<T>::largestBubbleID(bubMap.at(child.second).prevBubIDs, bubMap);
                    if (largestParent == curBubID)
                    {
                        lastBub = false;
                        curBubID = child.second;
                        curBub = &bubMap.at(curBubID);
                        break;
                    }
                }
            }
            if (lastBub) { break; }
        }
        // End this Oscillator; next, we decide whether or not to keep it
        osc.endTime = curBub->endTime;

        // Filter out if there is not enough solve data
        if (solveTimes.size() < 1) { continue; }
        // Filter out high freqency Oscillators
        if (*std::max_element(s_w0.begin(), s_w0.end()) > 2 * M_PI * 18000.) { continue; }
        // Filter out short blips
        if (osc.endTime - osc.startTime < 3 * 2 * M_PI / s_w0[0]) { continue; }

        if (transientPeriods > 0. && transientGain < 1.)
        {
            double periodCount = (osc.endTime - osc.startTime) / (2 * M_PI / s_w0[0]);
            if (periodCount < transientPeriods)
            {
                for (T& weight : f_weights) { weight *= transientGain; }
                transientOscillators++;
            }
        }

        // Clip each impulse's cutoff to the time available before the next event
        // (next chain link, or this Oscillator's endTime for the last impulse).
        // Otherwise a short-lived bubble (e.g. a SPLIT child that MERGEs back a
        // few samples later, or any non-largest chain leaf shorter than r/6 ms)
        // sees its F(t) = env(t/cutoff) * w * t^2 abruptly zeroed at the
        // coupled->uncoupled transition (Integrators.cpp _computeKCF only writes
        // _Fvals for the first _N_coupled oscillators) or at the chain switch in
        // _computeKCF when time crosses _forceData2(0, i). With the smoothstep
        // envelope, ending exactly at t=cutoff is C^1-zero, so this guarantees a
        // click-free release and a click-free switch to the next impulse.
        for (size_t i = 0; i < forceTimes.size(); i++)
        {
            const double next = (i + 1 < forceTimes.size())
                ? forceTimes[i + 1]
                : osc.endTime;
            const T max_dur = T(next - forceTimes[i]);
            if (max_dur > T(0) && f_cutoffs[i] > max_dur)
            {
                f_cutoffs[i] = max_dur;
            }
        }


        // Transfer solve data from temporary buffers to this Oscillator
        for (int i = 0; i < solveTimes.size(); i++)
        {
            Cvals.push_back(2. * Oscillator<T>::calcBeta(s_radii[i], s_w0[i], T(dampingCoeff)));
        }
        osc.solveTimes = solveTimes;
        osc.solveData.resize(6, solveTimes.size());

        osc.solveData.row(0) = Eigen::Map<Eigen::VectorX<T>>(s_radii.data(), s_radii.size());
        osc.solveData.row(1) = Eigen::Map<Eigen::VectorX<T>>(s_w0.data(), s_w0.size());
        osc.solveData.row(2) = Eigen::Map<Eigen::VectorX<T>>(s_x.data(), s_x.size());
        osc.solveData.row(3) = Eigen::Map<Eigen::VectorX<T>>(s_y.data(), s_y.size());
        osc.solveData.row(4) = Eigen::Map<Eigen::VectorX<T>>(s_z.data(), s_z.size());
        osc.solveData.row(5) = Eigen::Map<Eigen::VectorX<T>>(Cvals.data(), Cvals.size());

        osc.forceData.resize(3, forceTimes.size());
        osc.forceData.row(0) = Eigen::Map<Eigen::VectorX<T>>(forceTimes.data(), forceTimes.size());
        osc.forceData.row(1) = Eigen::Map<Eigen::VectorX<T>>(f_cutoffs.data(), f_cutoffs.size());
        osc.forceData.row(2) = Eigen::Map<Eigen::VectorX<T>>(f_weights.data(), f_weights.size());


        // Finally, add this Oscillator
        _oscillators.push_back(osc);
    } // END for (const std::pair<int, Bubble>& bubPair : bubMap)

    std::sort(_oscillators.begin(), _oscillators.end());

    if (timeJitterHalfWidth > 0.)
    {
        std::mt19937 gen;
        if (timeJitterSeed != 0ULL)
            gen.seed(static_cast<std::mt19937::result_type>(timeJitterSeed));
        else
        {
            std::random_device rd;
            gen.seed(rd());
        }
        std::uniform_real_distribution<double> dist(-timeJitterHalfWidth, timeJitterHalfWidth);
        for (auto& osc : _oscillators)
        {
            const double d = dist(gen);
            osc.startTime += d;
            osc.endTime += d;
            for (double& t : osc.solveTimes) { t += d; }
            for (int c = 0; c < osc.forceData.cols(); ++c)
                osc.forceData(0, c) += static_cast<T>(d);
        }
        std::sort(_oscillators.begin(), _oscillators.end());
    }

    std::set<double> eventTimesSet;
    for (const auto& osc : _oscillators)
    {
        eventTimesSet.insert(osc.startTime);
        eventTimesSet.insert(osc.endTime);
        if (denseEvents)
        {
            // Insert every per-sample-line solve time (after any jitter shift) so the
            //  integrator's K=w0^2 ramp follows the trackedBubInfo frequency column
            //  faithfully instead of linearly interpolating between only the first and
            //  last solveData column over the entire oscillator lifetime.
            for (double t : osc.solveTimes)
            {
                if (t > osc.startTime && t < osc.endTime)
                    eventTimesSet.insert(t);
            }
        }
    }
    _eventTimes.assign(eventTimesSet.begin(), eventTimesSet.end());

    if (transientPeriods > 0. && transientGain < 1.)
    {
        std::cout << "Transient oscillators attenuated = " << transientOscillators
                  << " (duration < " << transientPeriods
                  << " periods, gain=" << transientGain << ")" << std::endl;
    }

    // Damping (beta) summary across all surviving oscillators. We always print
    //  this so users can see the actual integrator damping, and so the effect
    //  of --damping-coeff is visible at runtime. solveData row 5 stores 2*beta.
    if (!_oscillators.empty())
    {
        double bMin = std::numeric_limits<double>::infinity();
        double bMax = -std::numeric_limits<double>::infinity();
        double bSum = 0.;
        long long nB = 0;
        for (const auto& osc : _oscillators)
        {
            const Eigen::Index N = osc.solveData.cols();
            for (Eigen::Index i = 0; i < N; ++i)
            {
                const double beta = 0.5 * static_cast<double>(osc.solveData(5, i));
                if (beta < bMin) bMin = beta;
                if (beta > bMax) bMax = beta;
                bSum += beta;
                ++nB;
            }
        }
        const double bMean = (nB > 0) ? (bSum / static_cast<double>(nB)) : 0.0;
        const double tauMin = (bMax > 0) ? (1.0 / bMax) : std::numeric_limits<double>::infinity();
        const double tauMax = (bMin > 0) ? (1.0 / bMin) : std::numeric_limits<double>::infinity();
        std::cout << "Damping coefficient applied = " << dampingCoeff
                  << "  (1.0 = unmodified Czerski/Deane radiative+viscous+thermal)" << std::endl;
        std::cout << "  beta range across all " << nB << " sample lines: ["
                  << bMin << ", " << bMax << "] s^-1, mean = " << bMean << " s^-1" << std::endl;
        std::cout << "  ringdown time tau = 1/beta in [" << tauMin << ", " << tauMax << "] s" << std::endl;

        if (_oscillators.size() == 1)
        {
            const auto& osc = _oscillators[0];
            const Eigen::Index N = osc.solveData.cols();
            const double w0First = static_cast<double>(osc.solveData(1, 0));
            const double w0Last = static_cast<double>(osc.solveData(1, N - 1));
            const double bFirst = 0.5 * static_cast<double>(osc.solveData(5, 0));
            const double bLast = 0.5 * static_cast<double>(osc.solveData(5, N - 1));
            std::cout << "  single oscillator: f_first=" << w0First / (2. * M_PI)
                      << " Hz beta_first=" << bFirst << " s^-1; f_last="
                      << w0Last / (2. * M_PI) << " Hz beta_last=" << bLast
                      << " s^-1" << std::endl;
        }
    }
}

template class Solver<float>;
template class Solver<double>;

} // namespace FluidSound