#include <chrono>
#include "sdsl/memory_management.hpp"

using namespace std::chrono;

namespace sdsl
{

template<>
void write_mem_log<CSV>(std::ostream& out,const memory_monitor& m)
{
    // write header
    out << "timestamp;memory_usage;event" << std::endl;

    auto first_ts = m.mem_events[0].timestamp;
    auto cur_event = m.events[0];
for (const auto& event : m.mem_events) {
        out << duration_cast<milliseconds>(event.timestamp-first_ts).count() << ";" << event.usage << ";";
    }
}

template<>
void write_mem_log<JSON>(std::ostream& out,const memory_monitor& m)
{
    // write header
    out << "[";

    auto first_ts = m.mem_events[0].timestamp;
    auto cur_event = m.events[0];
    for (size_t i=0; i<m.mem_events.size()-1; i++) {
        const auto& event = m.mem_events[i];
        out << "[" << duration_cast<milliseconds>(event.timestamp-first_ts).count() << "," << event.usage << "], ";
    }
    const auto& event = m.mem_events[m.mem_events.size()-1];
    out << "[" << duration_cast<milliseconds>(event.timestamp-first_ts).count() << "," << event.usage << "] ";
    out << "]";
}

}
