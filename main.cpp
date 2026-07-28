/** (c) 2024 Kangrui Xue
 *
 * \file main.cpp
 */

#include "FluidSound.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

typedef double precision;

/** */
void Run(const std::string& bubFile, const std::string& filteredFile, const std::string& outputFile, int srate, int scheme,
    double maxTime = 0., double timeJitterHalfWidth = 0., unsigned long long timeJitterSeed = 0ULL,
    double transientPeriods = 0., double transientGain = 1., double forcingCutoff = 0.0006,
    FluidSound::ForcingEnvelope forcingEnvelope = FluidSound::ForcingEnvelope::SMOOTHSTEP,
    bool denseEvents = false, double dampingCoeff = 1.0,
    bool applyListenerAttenuation = false,
    double listenerX = 0., double listenerY = 0., double listenerZ = 0., double listenerEpsilon = 1e-6,
    const std::string& energyLogFile = std::string(), const std::string& eventLogFile = std::string(),
    int energyStride = 48)
{
    std::cout << "runFluidSound: bubFile=\"" << bubFile << "\"";
    if (!filteredFile.empty())
        std::cout << " filteredFile=\"" << filteredFile << "\"";
    std::cout << " output=\"" << outputFile << "\" samplerate=" << srate << " scheme=" << scheme;
    if (maxTime > 0.)
        std::cout << " max_time=" << maxTime;
    if (timeJitterHalfWidth > 0.)
        std::cout << " time_jitter=+/- " << timeJitterHalfWidth << " s";
    if (transientPeriods > 0. && transientGain < 1.)
        std::cout << " transient_periods=" << transientPeriods << " transient_gain=" << transientGain;
    std::cout << " forcing_cutoff=" << forcingCutoff;
    std::cout << " forcing_envelope=" << (forcingEnvelope == FluidSound::ForcingEnvelope::SMOOTHSTEP ? "smoothstep" : "hard");
    std::cout << " dense_events=" << (denseEvents ? "on" : "off");
    std::cout << " damping_coeff=" << dampingCoeff;
    if (applyListenerAttenuation)
    {
        std::cout << " listener_pos=(" << listenerX << ", " << listenerY << ", " << listenerZ << ")"
                  << " listener_eps=" << listenerEpsilon;
    }
    std::cout << std::endl;

    // Check input file exists
    std::ifstream check(bubFile);
    if (!check.good())
    {
        std::cerr << "Error: Cannot open bubble file \"" << bubFile << "\": " << std::strerror(errno) << std::endl;
        throw std::runtime_error("Bubble file not found or not readable");
    }
    check.close();

    if (!filteredFile.empty())
    {
        std::ifstream fcheck(filteredFile);
        if (!fcheck.good())
        {
            std::cerr << "Error: Cannot open filtered file \"" << filteredFile << "\": " << std::strerror(errno) << std::endl;
            throw std::runtime_error("Filtered file not found or not readable");
        }
        fcheck.close();
    }

    std::cout << "Creating solver..." << std::endl;
    double dt = 1. / srate;
    FluidSound::Solver<precision> solver(bubFile, filteredFile, dt, scheme, 0., timeJitterHalfWidth, timeJitterSeed,
        transientPeriods, transientGain, forcingCutoff, forcingEnvelope, denseEvents, dampingCoeff,
        applyListenerAttenuation, listenerX, listenerY, listenerZ, listenerEpsilon);

    if (!energyLogFile.empty() || !eventLogFile.empty())
    {
        // If only one of the two log paths was given, derive the other next to it.
        std::string energyPath = energyLogFile.empty() ? (eventLogFile + "_energy.csv") : energyLogFile;
        std::string eventPath = eventLogFile.empty() ? (energyLogFile + "_events.csv") : eventLogFile;
        solver.enableEnergyLogging(energyPath, eventPath, energyStride);
    }

    const auto& evTimes = solver.eventTimes();
    if (evTimes.empty())
    {
        std::cerr << "Error: No event times (no bubbles loaded or all filtered out). Check bubble file format." << std::endl;
        throw std::runtime_error("Empty event times");
    }
    std::cout << "Simulation time range: " << evTimes.front() << " s to " << evTimes.back() << " s" << std::endl;

    double t_end = evTimes.back();
    if (maxTime > 0.)
    {
        t_end = std::min(t_end, maxTime);
        std::cout << "Time integration end: " << t_end << " s";
        if (t_end < evTimes.back())
            std::cout << " (early stop; bubble timeline ends at " << evTimes.back() << " s)";
        std::cout << std::endl;
    }

    std::ofstream out_file(outputFile);
    if (!out_file.good())
    {
        std::cerr << "Error: Cannot open \"" << outputFile << "\" for writing: " << std::strerror(errno) << std::endl;
        throw std::runtime_error("Cannot create output file");
    }

    // Simulate from t = 0. to t_end, logging sum of oscillators
    int tdx = 0;
    for (double t = 0.; t <= t_end; t += dt, tdx++)
    {
        out_file << solver.step() << std::endl;
        if (tdx % 9600 == 0) std::cout << "At time t = " << t << std::endl;
    }
    solver.printTimings();
    std::cout << "Done. Output written to " << outputFile << std::endl;
}

/** */
int main(int argc, char* argv[])
{
    try
    {
        // Usage: runFluidSound <bubFile> [filteredFile] <samplerate> <scheme> [-o output.txt] [--max-time SEC]
        //   bubFile:      full bubble graph (required for merge/split references)
        //   filteredFile: optional; if given, only these bubbles contribute to output
        //   samplerate:   48000 typical
        //   scheme:       0=uncoupled, 1=coupled
        //   -o/--output:  output waveform file (default: output.txt)
        //   --max-time:   if > 0, stop simulation at this time (seconds) and write partial waveform
        //   --time-jitter W: each oscillator shifted by U(-W,W) seconds (desync grid clicks)
        //   --time-jitter-seed: optional RNG seed (default: random_device)
        //   --transient-periods N: attenuate start forcing for oscillators lasting fewer than N periods
        //   --transient-gain G: forcing multiplier for transient oscillators (default 1 = unchanged)
        //   --forcing-cutoff SEC: max forcing duration for start impulses (default 0.0006)
        //   --forcing-envelope hard|smoothstep: start impulse envelope (default hard)
        //   --dense-events: insert each per-sample-line solveTime into the integrator
        //                   event-time set so K=w0^2 is sampled at every trackedBubInfo
        //                   row (instead of a linear ramp between only the first and
        //                   last solve column over each oscillator's lifetime). Default off.
        //   --damping-coeff C: multiplier on the per-sample beta computed by
        //                   Oscillator::calcBeta. Default 1.0 (Czerski/Deane). Use C<1
        //                   to lengthen ringdown (longer audible ring), C>1 to shorten it.
        //   --listener-position X Y Z: enable per-oscillator 1/distance-to-listener
        //                   attenuation on the audio sum. The listener sits at (X, Y, Z)
        //                   in the same coordinate space as the bubble positions in the
        //                   input trackedBubInfo. The factor 1 / max(d, eps) is applied
        //                   AFTER integration (coupled dynamics untouched) per oscillator
        //                   at every audio sample. eps defaults to 1e-6 and can be
        //                   overridden via --listener-epsilon. Default off (no attenuation).
        //   --listener-epsilon E: lower clamp on d used by --listener-position to avoid
        //                   division blow-ups when a bubble sits exactly on the listener.
        //   --energy-log F: write population oscillator-energy CSV (time,E_tot,KE,PE,
        //                   P_in,P_diss,cumulative work by event type, ...) to F.
        //                   Audit-only; the waveform is unchanged.
        //   --event-log F:  write one CSV row per forcing impulse (event type, radius,
        //                   w0, cutoff, weight, E_before, E_after) to F.
        //   --energy-stride N: energy CSV row every N samples (default 48 = 1 kHz @48k).
        std::string bubFile("../Scenes/GlassPour/trackedBubInfo.txt");
        std::string filteredFile;
        std::string outputFile("output.txt");
        int samplerate = 48000;
        int scheme = 1;
        double maxTime = 0.;
        double timeJitterHalfWidth = 0.;
        unsigned long long timeJitterSeed = 0ULL;
        double transientPeriods = 0.;
        double transientGain = 1.;
        double forcingCutoff = 0.0006;
        FluidSound::ForcingEnvelope forcingEnvelope = FluidSound::ForcingEnvelope::SMOOTHSTEP;
        bool denseEvents = false;
        double dampingCoeff = 1.0;
        bool applyListenerAttenuation = false;
        double listenerX = 0., listenerY = 0., listenerZ = 0.;
        double listenerEpsilon = 1e-6;
        std::string energyLogFile;
        std::string eventLogFile;
        int energyStride = 48;

        std::vector<std::string> posArgs;
        for (int i = 1; i < argc; i++)
        {
            std::string arg(argv[i]);
            if (arg == "-o" || arg == "--output")
            {
                if (i + 1 < argc) { outputFile = std::string(argv[++i]); }
            }
            else if (arg == "--max-time" || arg == "--max_time")
            {
                if (i + 1 < argc) { maxTime = std::atof(argv[++i]); }
            }
            else if (arg == "--time-jitter" || arg == "--time_jitter")
            {
                if (i + 1 < argc) { timeJitterHalfWidth = std::atof(argv[++i]); }
            }
            else if (arg == "--time-jitter-seed" || arg == "--time_jitter_seed")
            {
                if (i + 1 < argc) { timeJitterSeed = std::strtoull(argv[++i], nullptr, 10); }
            }
            else if (arg == "--transient-periods" || arg == "--transient_periods")
            {
                if (i + 1 < argc) { transientPeriods = std::atof(argv[++i]); }
            }
            else if (arg == "--transient-gain" || arg == "--transient_gain")
            {
                if (i + 1 < argc) { transientGain = std::atof(argv[++i]); }
            }
            else if (arg == "--forcing-cutoff" || arg == "--forcing_cutoff")
            {
                if (i + 1 < argc) { forcingCutoff = std::atof(argv[++i]); }
            }
            else if (arg == "--forcing-envelope" || arg == "--forcing_envelope")
            {
                if (i + 1 < argc)
                {
                    std::string envelope = std::string(argv[++i]);
                    if (envelope == "smoothstep") { forcingEnvelope = FluidSound::ForcingEnvelope::SMOOTHSTEP; }
                    else if (envelope == "hard") { forcingEnvelope = FluidSound::ForcingEnvelope::HARD; }
                    else { throw std::runtime_error("Invalid --forcing-envelope (expected hard or smoothstep): " + envelope); }
                }
            }
            else if (arg == "--dense-events" || arg == "--dense_events")
            {
                denseEvents = true;
            }
            else if (arg == "--no-dense-events" || arg == "--no_dense_events")
            {
                denseEvents = false;
            }
            else if (arg == "--damping-coeff" || arg == "--damping_coeff")
            {
                if (i + 1 < argc) { dampingCoeff = std::atof(argv[++i]); }
            }
            else if (arg == "--listener-position" || arg == "--listener_position")
            {
                if (i + 3 < argc)
                {
                    listenerX = std::atof(argv[++i]);
                    listenerY = std::atof(argv[++i]);
                    listenerZ = std::atof(argv[++i]);
                    applyListenerAttenuation = true;
                }
                else
                {
                    throw std::runtime_error("--listener-position requires 3 arguments: X Y Z");
                }
            }
            else if (arg == "--listener-epsilon" || arg == "--listener_epsilon")
            {
                if (i + 1 < argc) { listenerEpsilon = std::atof(argv[++i]); }
            }
            else if (arg == "--energy-log" || arg == "--energy_log")
            {
                if (i + 1 < argc) { energyLogFile = std::string(argv[++i]); }
            }
            else if (arg == "--event-log" || arg == "--event_log")
            {
                if (i + 1 < argc) { eventLogFile = std::string(argv[++i]); }
            }
            else if (arg == "--energy-stride" || arg == "--energy_stride")
            {
                if (i + 1 < argc) { energyStride = std::atoi(argv[++i]); }
            }
            else
            {
                posArgs.push_back(arg);
            }
        }

        if (posArgs.size() > 0) bubFile = posArgs[0];
        if (posArgs.size() >= 4)
        {
            filteredFile = posArgs[1];
            samplerate = std::atoi(posArgs[2].c_str());
            scheme = std::atoi(posArgs[3].c_str());
        }
        else if (posArgs.size() >= 3)
        {
            samplerate = std::atoi(posArgs[1].c_str());
            scheme = std::atoi(posArgs[2].c_str());
        }
        else if (posArgs.size() >= 2)
        {
            samplerate = std::atoi(posArgs[1].c_str());
        }

        Run(bubFile, filteredFile, outputFile, samplerate, scheme, maxTime, timeJitterHalfWidth, timeJitterSeed,
            transientPeriods, transientGain, forcingCutoff, forcingEnvelope, denseEvents, dampingCoeff,
            applyListenerAttenuation, listenerX, listenerY, listenerZ, listenerEpsilon,
            energyLogFile, eventLogFile, energyStride);
        return 0;
    }
    catch (const std::out_of_range& e)
    {
        std::cerr << "Error (missing bubble reference): " << e.what() << std::endl;
        std::cerr << "Hint: Filtered bubble files may have broken merge/split chains. "
                  << "Each bubble's prevBubIDs and nextBubIDs must exist in the file." << std::endl;
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Error: Unknown exception" << std::endl;
        return 1;
    }
}