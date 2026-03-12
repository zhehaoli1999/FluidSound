/** (c) 2024 Kangrui Xue
 *
 * \file BubbleUtils.cpp
 * \brief Implements BubbleUtils class
 */

#include "BubbleUtils.h"


namespace FluidSound {

/** */
template <typename T>
void BubbleUtils<T>::_parseBubble(std::pair<int, Bubble<T>> &bubPair, std::ifstream& in)
{
    Bubble<T>& bub = bubPair.second;
    std::string line; std::istringstream is;

    int bubID; double time;
    T freqHz, x, y, z;
    T pressure = 101450.;  // if pressure not specified, assume default value of 101450 Pa

    // Line 1: 'Bub <unique bubble ID> <radius>'
    std::getline(in, line);
    if (line.empty())
    {
        bubPair.first = -1;
        return;
    }
    is = std::istringstream(line.substr(4));
    is >> bub.bubID; bubPair.first = bub.bubID;
    is >> bub.radius;

    // Line 2: '  Start: <event type> <start time> <previous bubble ID(s)>'
    std::getline(in, line);
    if (line.size() < 10)
    {
        throw std::runtime_error("Invalid Start line (too short): " + line);
    }
    switch (line[9])
    {
        case 'N': bub.startType = EventType::ENTRAIN; break;
        case 'M': bub.startType = EventType::MERGE; break;
        case 'S': bub.startType = EventType::SPLIT; break;
        default: throw std::runtime_error(std::string("Invalid bubble start event type (expected N/M/S), got '") + line[9] + "': " + line);
    }
    is = std::istringstream(line.substr(11));
    is >> bub.startTime;
    while (is >> bubID)
    {
        bub.prevBubIDs.push_back(bubID);
    }

    // Line 3-n: '  <time> <freqHz> <x> <y> <z> <pressure inside bubble (optional)>'
    std::getline(in, line);
    while (line[2] != 'E' && line[2] != 'B' && !in.eof())
    {
        pressure = 101450.;  // default if not in line
        is = std::istringstream(line);
        is >> time >> freqHz >> x >> y >> z >> pressure;
        std::getline(in, line);

        bub.solveTimes.push_back(time); 
        bub.w0.push_back(2 * M_PI * freqHz);
        bub.x.push_back(x); bub.y.push_back(y); bub.z.push_back(z);
        bub.pressure.push_back(pressure);
    };

    // Line n: '  End: <event type> <end time> <next bubble ID(s)>'...
    if (line[2] == 'E')
    {
        if (line.size() < 8)
        {
            throw std::runtime_error("Invalid End line (too short): " + line);
        }
        switch (line[7])
        {
            case 'M': bub.endType = EventType::MERGE; break;
            case 'S': bub.endType = EventType::SPLIT; break;
            case 'C': bub.endType = EventType::COLLAPSE; break;
            default: throw std::runtime_error(std::string("Invalid bubble end event type (expected M/S/C), got '") + line[7] + "': " + line);
        }
        is = std::istringstream(line.substr(9));
        is >> bub.endTime;
        while (is >> bubID)
        {
            bub.nextBubIDs.push_back(bubID);
        }

        if (bub.endType == EventType::MERGE && bub.nextBubIDs.size() != 1)
        {
            throw std::runtime_error("Bubble " + std::to_string(bub.bubID) + " did not merge to one bubble (got " + std::to_string(bub.nextBubIDs.size()) + " next IDs)");
        }
        else if (bub.endType == EventType::SPLIT && bub.nextBubIDs.size() < 2)
        {
            throw std::runtime_error("Bubble " + std::to_string(bub.bubID) + " split to less than 2 bubbles (got " + std::to_string(bub.nextBubIDs.size()) + " next IDs)");
        }
    }
    // ...however, if end info not present, use default values
    else if (line[2] == 'B') 
    {
        if (bub.solveTimes.empty()) { bub.endTime = bub.startTime + 0.001; }
        else { bub.endTime = bub.solveTimes.back(); }
        bub.endType = EventType::COLLAPSE;
    }
}

/** */
template <typename T>
void BubbleUtils<T>::loadBubbleFile(std::map<int, Bubble<T>>& bubMap, const std::string& bubFile)
{
    bubMap.clear();
    std::ifstream inFile(bubFile);
    if (!inFile.is_open() || !inFile.good())
    {
        throw std::runtime_error(std::string("Failed to open bubble file: ") + bubFile);
    }
    while (inFile.good())
    {
        std::pair<int, Bubble<T>> bubPair;
        _parseBubble(bubPair, inFile);
        if (bubPair.first != -1) { bubMap.insert(bubPair); }
    }
}

/** */
template <typename T>
int BubbleUtils<T>::largestBubbleID(const std::vector<int>& bubIDs, const std::map<int, Bubble<T>>& bubMap)
{
    T max_radius = 0.;
    int max_bubID = 0;
    for (int id : bubIDs)
    {
        if (bubMap.at(id).radius > max_radius)
        {
            max_radius = bubMap.at(id).radius;
            max_bubID = id;
        }
    }
    return max_bubID;
}

/** */
std::set<int> parseBubbleIDsFromFile(const std::string& bubFile)
{
    std::set<int> ids;
    std::ifstream in(bubFile);
    if (!in.is_open() || !in.good())
    {
        throw std::runtime_error(std::string("Failed to open bubble ID file: ") + bubFile);
    }
    std::string line;
    while (std::getline(in, line))
    {
        if (line.size() >= 4 && line.substr(0, 4) == "Bub ")
        {
            int id;
            std::istringstream is(line.substr(4));
            if (is >> id)
            {
                ids.insert(id);
            }
        }
    }
    return ids;
}

template class BubbleUtils<float>;
template class BubbleUtils<double>;

} // namespace FluidSound