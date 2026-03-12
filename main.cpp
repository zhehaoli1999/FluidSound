/** (c) 2024 Kangrui Xue
 *
 * \file main.cpp
 */

#include "FluidSound.h"

#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

typedef double precision;

/** */
void Run(const std::string& bubFile, const std::string& filteredFile, const std::string& outputFile, int srate, int scheme)
{
    std::cout << "runFluidSound: bubFile=\"" << bubFile << "\"";
    if (!filteredFile.empty())
        std::cout << " filteredFile=\"" << filteredFile << "\"";
    std::cout << " output=\"" << outputFile << "\" samplerate=" << srate << " scheme=" << scheme << std::endl;

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
    FluidSound::Solver<precision> solver(bubFile, filteredFile, dt, scheme);

    const auto& evTimes = solver.eventTimes();
    if (evTimes.empty())
    {
        std::cerr << "Error: No event times (no bubbles loaded or all filtered out). Check bubble file format." << std::endl;
        throw std::runtime_error("Empty event times");
    }
    std::cout << "Simulation time range: " << evTimes.front() << " s to " << evTimes.back() << " s" << std::endl;

    std::ofstream out_file(outputFile);
    if (!out_file.good())
    {
        std::cerr << "Error: Cannot open \"" << outputFile << "\" for writing: " << std::strerror(errno) << std::endl;
        throw std::runtime_error("Cannot create output file");
    }

    // Simulate from t = 0. to t = 'tf' seconds, logging sum of oscillators
    int tdx = 0;
    for (double t = 0.; t <= evTimes.back(); t += dt, tdx++)
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
        // Usage: runFluidSound <bubFile> [filteredFile] <samplerate> <scheme> [-o output.txt]
        //   bubFile:      full bubble graph (required for merge/split references)
        //   filteredFile: optional; if given, only these bubbles contribute to output
        //   samplerate:   48000 typical
        //   scheme:       0=uncoupled, 1=coupled
        //   -o/--output:  output waveform file (default: output.txt)
        std::string bubFile("../Scenes/GlassPour/trackedBubInfo.txt");
        std::string filteredFile;
        std::string outputFile("output.txt");
        int samplerate = 48000;
        int scheme = 1;

        std::vector<std::string> posArgs;
        for (int i = 1; i < argc; i++)
        {
            std::string arg(argv[i]);
            if (arg == "-o" || arg == "--output")
            {
                if (i + 1 < argc) { outputFile = std::string(argv[++i]); }
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

        Run(bubFile, filteredFile, outputFile, samplerate, scheme);
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