/**
 * Create trackedBubInfo variants with frozen or Minnaert frequencies.
 *
 * No amplitude/damping computation. Only modifies frequencies and writes
 * original trackedBubInfo format (time freq x y z pressure) for runFluidSound.
 *
 * Usage:
 *   create_freq_variants <input.txt> --frozen-freq <out.txt> [--minnaert-freq <out.txt>]
 *   create_freq_variants trackedBubInfo.txt --frozen-freq trackedBubInfo_frozen_freq.txt --minnaert-freq trackedBubInfo_minnaert_freq.txt
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#include "BubbleUtils.h"

/** Minnaert frequency: f = 3.26 / r (Hz), r in meters */
static const double MINNAERT_COEFF = 3.283243423687599;

/** Write trackedBubInfo in original format (time freq x y z pressure) with modified frequencies.
 *  freqMode: "frozen" = first frequency per bubble, "minnaert" = 3.26/radius per bubble */
static void writeFreqVariant(const std::string& outPath,
    const std::map<int, FluidSound::Bubble<double>>& bubMap,
    const std::string& freqMode)
{
    std::ofstream out(outPath);
    if (!out) throw std::runtime_error("Cannot open output file: " + outPath);

    std::vector<int> sortedIDs;
    for (const auto& pair : bubMap) sortedIDs.push_back(pair.first);
    std::sort(sortedIDs.begin(), sortedIDs.end());  // same order as original: by bubble ID

    out << std::setprecision(12);

    for (int bubID : sortedIDs)
    {
        const FluidSound::Bubble<double>& bub = bubMap.at(bubID);
        double frozenFreqHz;
        if (freqMode == "frozen")
        {
            if (bub.w0.empty())
                frozenFreqHz = MINNAERT_COEFF / bub.radius;
            else
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
            case FluidSound::EventType::SPLIT:  out << "S"; break;
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
        else if (inFile == "trackedBubInfo.txt")
        {
            inFile = arg;
        }
    }

    if (frozenFreqOut.empty() && minnaertFreqOut.empty())
    {
        std::cerr << "Usage: create_freq_variants <input.txt> --frozen-freq <out.txt> [--minnaert-freq <out.txt>]\n";
        return 1;
    }

    std::cout << "Reading " << inFile << "..." << std::endl;
    std::map<int, FluidSound::Bubble<double>> bubMap;
    FluidSound::BubbleUtils<double>::loadBubbleFile(bubMap, inFile);
    std::cout << "Loaded " << bubMap.size() << " bubbles." << std::endl;

    if (!frozenFreqOut.empty())
    {
        std::cout << "Writing frozen-frequency variant to " << frozenFreqOut << "..." << std::endl;
        writeFreqVariant(frozenFreqOut, bubMap, "frozen");
        std::cout << "Done (frozen freq)." << std::endl;
    }
    if (!minnaertFreqOut.empty())
    {
        std::cout << "Writing Minnaert-frequency variant to " << minnaertFreqOut << "..." << std::endl;
        writeFreqVariant(minnaertFreqOut, bubMap, "minnaert");
        std::cout << "Done (Minnaert freq)." << std::endl;
    }

    std::cout << "Done." << std::endl;
    return 0;
}
