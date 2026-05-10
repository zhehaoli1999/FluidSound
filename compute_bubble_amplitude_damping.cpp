/**
 * Standalone tool to compute amplitude and damping for each bubble sample
 * from trackedBubInfo.txt, and write an updated file with these fields.
 *
 * Uses FluidSound's Oscillator::calcBeta for damping (2β) and forcing models
 * (CzerskiJetForcing, MergeForcing) for amplitude, following [Langlois et al. 2016].
 *
 * Output format: same as input, with weight, cutoff, and damping appended per data line:
 *   time freqHz x y z [pressure] weight cutoff damping
 *
 * Optional modes for frequency experiment (full vs frozen vs Minnaert):
 *   --frozen-freq <path>   Write trackedBubInfo with first (initial) frequency per bubble
 *   --minnaert-freq <path> Write trackedBubInfo with Minnaert frequency (3.26/r Hz) per bubble
 *   These output original format (time freq x y z pressure) for runFluidSound.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

#include "BubbleUtils.h"
#include "Oscillator.h"

// Deterministic MergeForcing (fixed frac for reproducibility; FluidSound uses random)
template <typename T>
static std::pair<T, T> mergeForcingDeterministic(T radius, T r1, T r2, T frac = 0.6)
{
    const T RHO_WATER = 998.;
    const T GAMMA = 1.4;
    const T SIGMA = 0.0726;
    const T ATM = 101325;
    const T MAX_CUTOFF = 0.0006;

    T factor = std::pow(2. * SIGMA * r1 * r2 / (RHO_WATER * (r1 + r2)), 0.25);
    T cutoff = std::min(MAX_CUTOFF, 0.5 / (3. / radius));
    T tmp = std::pow(frac * std::min(r1, r2) / 2. / factor, 2);
    cutoff = std::min(cutoff, tmp);

    T pressure_in0 = (ATM + 2. * SIGMA / radius);
    T weight = 6. * SIGMA * GAMMA * pressure_in0 / (RHO_WATER * radius * radius * radius);
    T mass = (RHO_WATER / (4. * M_PI * radius));

    return std::pair<T, T>(cutoff, weight / mass);
}

// Compute forcing (weight, cutoff) for a bubble based on its start event
template <typename T>
std::pair<T, T> computeForcing(const FluidSound::Bubble<T>& bub, const std::map<int, FluidSound::Bubble<T>>& bubMap)
{
    using namespace FluidSound;

    if (bub.startType == EventType::ENTRAIN)
    {
        auto f = Oscillator<T>::CzerskiJetForcing(bub.radius);
        return std::make_pair(f.second, f.first);
    }
    if (bub.startType == EventType::SPLIT && !bub.prevBubIDs.empty())
    {
        int parentID = bub.prevBubIDs[0];
        if (bubMap.count(parentID) && bubMap.at(parentID).radius >= bub.radius)
        {
            auto f = Oscillator<T>::CzerskiJetForcing(bub.radius);
            return std::make_pair(f.second, f.first);
        }
    }
    if (bub.startType == EventType::MERGE && bub.prevBubIDs.size() == 2)
    {
        int p1 = bub.prevBubIDs[0], p2 = bub.prevBubIDs[1];
        if (bubMap.count(p1) && bubMap.count(p2))
        {
            bool allMerge = (bubMap.at(p1).endType == EventType::MERGE &&
                            bubMap.at(p2).endType == EventType::MERGE);
            if (allMerge)
            {
                T r1 = bubMap.at(p1).radius, r2 = bubMap.at(p2).radius;
                if (r1 + r2 > bub.radius)
                {
                    T v1 = 4. / 3. * M_PI * r1 * r1 * r1;
                    T v2 = 4. / 3. * M_PI * r2 * r2 * r2;
                    T vn = 4. / 3. * M_PI * bub.radius * bub.radius * bub.radius;
                    T diff = v1 + v2 - vn;
                    if (diff <= std::max(v1, v2))
                    {
                        if (v1 > v2) v1 -= diff; else v2 -= diff;
                        r1 = std::pow(3. / 4. / M_PI * v1, 1. / 3.);
                        r2 = std::pow(3. / 4. / M_PI * v2, 1. / 3.);
                    }
                }
                auto f = mergeForcingDeterministic(bub.radius, r1, r2);
                return std::make_pair(f.second, f.first);
            }
        }
    }
    return std::make_pair(T(0.), T(0.));  // no forcing
}

static void processAndWrite(const std::string& outPath,
                            const std::map<int, FluidSound::Bubble<double>>& bubMap,
                            double tMin, double tMax)
{
    std::ofstream out(outPath);
    if (!out)
    {
        throw std::runtime_error("Cannot open output file");
    }

    // Sort by min solve time for time-ordered progress
    std::vector<int> sortedIDs;
    for (const auto& pair : bubMap) sortedIDs.push_back(pair.first);
    std::sort(sortedIDs.begin(), sortedIDs.end(), [&bubMap](int a, int b) {
        double ta = bubMap.at(a).solveTimes.empty() ? bubMap.at(a).startTime : bubMap.at(a).solveTimes.front();
        double tb = bubMap.at(b).solveTimes.empty() ? bubMap.at(b).startTime : bubMap.at(b).solveTimes.front();
        return ta < tb;
    });

    const size_t total = sortedIDs.size();
    double tCurrent = tMin;

    for (size_t idx = 0; idx < total; idx++)
    {
        int bubID = sortedIDs[idx];
        const FluidSound::Bubble<double>& bub = bubMap.at(bubID);
        auto forcing = computeForcing(bub, bubMap);
        double weight = forcing.first;
        double cutoff = forcing.second;

        out << "Bub " << bub.bubID << " " << bub.radius << "\n";
        out << "  Start: ";
        switch (bub.startType)
        {
            case FluidSound::EventType::ENTRAIN: out << "N"; break;
            case FluidSound::EventType::MERGE:   out << "M"; break;
            case FluidSound::EventType::SPLIT:   out << "S"; break;
        }
        out << " " << bub.startTime;
        for (int id : bub.prevBubIDs) out << " " << id;
        out << "\n";

        for (size_t i = 0; i < bub.solveTimes.size(); i++)
        {
            double w0 = bub.w0[i];
            double freqHz = w0 / (2. * M_PI);
            double damping = 2. * FluidSound::Oscillator<double>::calcBeta(bub.radius, w0);
            out << "  " << bub.solveTimes[i] << " " << freqHz << " " << bub.x[i] << " "
                << bub.y[i] << " " << bub.z[i] << " "
                << bub.pressure[i] << " "
                << weight << " " << cutoff << " " << damping << "\n";
            if (bub.solveTimes[i] > tCurrent) tCurrent = bub.solveTimes[i];
        }
        if (!bub.solveTimes.empty() && bub.solveTimes.back() > tCurrent)
            tCurrent = bub.solveTimes.back();

        double pct = (tMax > tMin) ? 100. * (tCurrent - tMin) / (tMax - tMin) : 100.;
        std::cerr << "\r  [" << std::fixed << std::setprecision(1) << pct << "%] t=" << tCurrent << "s / " << tMax << "s    " << std::flush;

        out << "  End: ";
        switch (bub.endType)
        {
            case FluidSound::EventType::MERGE:   out << "M"; break;
            case FluidSound::EventType::SPLIT:   out << "S"; break;
            case FluidSound::EventType::COLLAPSE: out << "C"; break;
        }
        out << " " << bub.endTime;
        for (int id : bub.nextBubIDs) out << " " << id;
        out << "\n";
    }
}

/** Minnaert frequency: f = 3.26 / r (Hz), r in meters */
static const double MINNAERT_COEFF = 3.283243423687599;

/** Write trackedBubInfo in original format (time freq x y z pressure) with modified frequencies.
 *  freqMode: "frozen" = first frequency per bubble, "minnaert" = 3.26/radius per bubble */
static void writeTrackedBubInfoFreqVariant(const std::string& outPath,
    const std::map<int, FluidSound::Bubble<double>>& bubMap,
    const std::string& freqMode)
{
    std::ofstream out(outPath);
    if (!out) throw std::runtime_error("Cannot open output file: " + outPath);

    std::vector<int> sortedIDs;
    for (const auto& pair : bubMap) sortedIDs.push_back(pair.first);
    std::sort(sortedIDs.begin(), sortedIDs.end(), [&bubMap](int a, int b) {
        double ta = bubMap.at(a).solveTimes.empty() ? bubMap.at(a).startTime : bubMap.at(a).solveTimes.front();
        double tb = bubMap.at(b).solveTimes.empty() ? bubMap.at(b).startTime : bubMap.at(b).solveTimes.front();
        return ta < tb;
    });

    for (int bubID : sortedIDs)
    {
        const FluidSound::Bubble<double>& bub = bubMap.at(bubID);
        double frozenFreqHz;
        if (freqMode == "frozen")
        {
            if (bub.w0.empty()) continue;
            frozenFreqHz = bub.w0[0] / (2. * M_PI);
        }
        else if (freqMode == "minnaert")
        {
            frozenFreqHz = MINNAERT_COEFF / bub.radius;
        }
        else
        {
            throw std::runtime_error("Unknown freqMode: " + freqMode);
        }

        out << "Bub " << bub.bubID << " " << bub.radius << "\n";
        out << "  Start: ";
        switch (bub.startType)
        {
            case FluidSound::EventType::ENTRAIN: out << "N"; break;
            case FluidSound::EventType::MERGE:   out << "M"; break;
            case FluidSound::EventType::SPLIT:   out << "S"; break;
        }
        out << " " << bub.startTime;
        for (int id : bub.prevBubIDs) out << " " << id;
        out << "\n";

        for (size_t i = 0; i < bub.solveTimes.size(); i++)
        {
            out << "  " << bub.solveTimes[i] << " " << frozenFreqHz << " " << bub.x[i] << " "
                << bub.y[i] << " " << bub.z[i] << " " << bub.pressure[i] << "\n";
        }

        out << "  End: ";
        switch (bub.endType)
        {
            case FluidSound::EventType::MERGE:   out << "M"; break;
            case FluidSound::EventType::SPLIT:   out << "S"; break;
            case FluidSound::EventType::COLLAPSE: out << "C"; break;
        }
        out << " " << bub.endTime;
        for (int id : bub.nextBubIDs) out << " " << id;
        out << "\n";
    }
}

int main(int argc, char* argv[])
{
    std::string inFile = "trackedBubInfo.txt";
    std::string outFile = "trackedBubInfo_with_amplitude_damping.txt";
    std::string frozenFreqOut;
    std::string minnaertFreqOut;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--frozen-freq" && i + 1 < argc)
        {
            frozenFreqOut = argv[++i];
        }
        else if (arg == "--minnaert-freq" && i + 1 < argc)
        {
            minnaertFreqOut = argv[++i];
        }
        else if (i == 1) inFile = arg;
        else if (i == 2) outFile = arg;
    }

    std::cout << "Reading " << inFile << "..." << std::endl;
    std::map<int, FluidSound::Bubble<double>> bubMap;
    FluidSound::BubbleUtils<double>::loadBubbleFile(bubMap, inFile);
    std::cout << "Loaded " << bubMap.size() << " bubbles." << std::endl;

    double tMin = 1e99, tMax = -1e99;
    for (const auto& pair : bubMap)
    {
        for (double t : pair.second.solveTimes)
        {
            if (t < tMin) tMin = t;
            if (t > tMax) tMax = t;
        }
    }
    if (tMin > tMax) { tMin = 0.; tMax = 1.; }
    std::cout << "Simulation time range: " << tMin << " - " << tMax << " s" << std::endl;

    if (!frozenFreqOut.empty())
    {
        std::cout << "Writing frozen-frequency variant to " << frozenFreqOut << "..." << std::endl;
        writeTrackedBubInfoFreqVariant(frozenFreqOut, bubMap, "frozen");
        std::cout << "Done (frozen freq)." << std::endl;
    }
    if (!minnaertFreqOut.empty())
    {
        std::cout << "Writing Minnaert-frequency variant to " << minnaertFreqOut << "..." << std::endl;
        writeTrackedBubInfoFreqVariant(minnaertFreqOut, bubMap, "minnaert");
        std::cout << "Done (Minnaert freq)." << std::endl;
    }

    std::cout << "Computing amplitude and damping, writing to " << outFile << "..." << std::endl;
    processAndWrite(outFile, bubMap, tMin, tMax);
    std::cerr << "\r  [100.0%] t=" << tMax << "s / " << tMax << "s    " << std::endl;
    std::cout << "Done." << std::endl;

    return 0;
}
