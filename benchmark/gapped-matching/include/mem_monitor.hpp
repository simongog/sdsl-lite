// mem monitor
// Copyright (c) 2014, Matthias Petri, All rights reserved.

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
#pragma once

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <stdexcept>

#ifndef MEM_MON_MEM_LIMIT_MB
#define MEM_MON_MEM_LIMIT_MB 32
#endif

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;
using std::vector;

struct mem_stat {
    steady_clock::time_point timestamp;
    uint64_t pid;
    uint64_t VmPeak;
    uint64_t VmSize;
    uint64_t VmHWM;
    uint64_t VmRSS;
    uint64_t VmData;
    uint64_t VmPTE;
    uint64_t event_id;
};

class mem_monitor
{
    private:
        steady_clock::time_point                                    m_start_timestamp;
        vector<mem_stat>                                                      m_stats;
        vector<std::string>                                             m_event_names;
        std::ofstream                                                            m_os;
        milliseconds                                                    m_granularity;
        std::condition_variable                                                  m_cv;
        std::mutex                                                            m_mutex;
        std::thread                                                          m_thread;
        std::atomic<bool>                               m_run = ATOMIC_VAR_INIT(true);
        std::atomic<uint64_t>                     m_cur_event_id = ATOMIC_VAR_INIT(0);
        uint64_t                             m_memory_limit_mb = MEM_MON_MEM_LIMIT_MB;
        bool                                                           m_write_header;

    private:
        void monitor() {
            while (m_run) {
                m_stats.emplace_back(get_current_stats());

                if (write_needed()) {
                    flush();
                }

                if (m_run) {
                    std::unique_lock<std::mutex> lk(m_mutex);
                    m_cv.wait_for(lk, m_granularity);
                }
            }
        }

        uint64_t extract_number(std::string& line, size_t start_pos,
                                bool extension = false) {
            auto num_end = line.find_first_of(' ', start_pos);
            if (num_end == std::string::npos) {
                num_end = line.size() - 1;
            }
            uint64_t num = std::strtoull(line.c_str() + start_pos, NULL, 10);

            if (extension) {
                if (line.back() == 'B') {
                    if (line[line.size() - 2] == 'k' || line[line.size() - 2] == 'K') {
                        num *= 1024;
                    }
                    if (line[line.size() - 2] == 'm' || line[line.size() - 2] == 'M') {
                        num *= 1024 * 1024;
                    }
                    if (line[line.size() - 2] == 'g' || line[line.size() - 2] == 'G') {
                        num *= 1024 * 1024 * 1024;
                    }
                } else {
                    throw std::invalid_argument("no extension found during line parsing");
                }
            }
            return num;
        }

        bool write_needed() {
            if (m_memory_limit_mb * 1024 * 1024 < m_stats.size() * sizeof(mem_stat))
                return true;
            else
                return false;
        }

        void flush() {
            if (m_write_header) {
                m_os << "time_ms;pid;VmPeak;VmSize;VmHWM;VmRSS;VmData;VmPTE;event\n";
                m_write_header = false;
            }

            for (const auto &s : m_stats) {
                m_os << duration_cast<milliseconds>(s.timestamp - m_start_timestamp).count() + 1 << ";" // round up!
                     << s.pid    << ";"
                     << s.VmPeak << ";"
                     << s.VmSize << ";"
                     << s.VmHWM  << ";"
                     << s.VmRSS  << ";"
                     << s.VmData << ";"
                     << s.VmPTE  << ";"
                     << '"' << m_event_names[s.event_id] << '"' << "\n";
            }
            m_stats.clear();
        }

    public:
        // delete the constructors and operators we don't want
        mem_monitor() = delete;
        mem_monitor(const mem_monitor&) = delete;
        mem_monitor(mem_monitor&&) = delete;
        mem_monitor& operator=(const mem_monitor&) = delete;
        mem_monitor& operator=(mem_monitor&&) = delete;

        mem_monitor(const std::string& file_name,
                    milliseconds granularity = milliseconds(50))
            : m_os(file_name), m_granularity(granularity), m_write_header(true) {
            // some init stuff
            m_event_names.push_back(""); // default empty event
            m_start_timestamp = steady_clock::now();

            if (!m_os.is_open()) { // output stream open?
                throw std::ios_base::failure("memory monitor output file could not be opened.");
            }

            // spawn the thread
            m_thread = std::thread(&mem_monitor::monitor, this);
        }

        ~mem_monitor() {
            m_run = false;
            m_cv.notify_one(); // notify the sleeping thread
            m_thread.join();
            flush();
            m_os.close();
        }

        void event(const std::string& ev) {
            m_event_names.push_back(ev);
            m_cur_event_id++;
        }


        mem_stat get_current_stats() {
            mem_stat stat;
            stat.timestamp = steady_clock::now();
            stat.event_id = m_cur_event_id;

            // read memory stats
            {
                std::ifstream pfs("/proc/self/status");
                std::string line;
                while (std::getline(pfs, line)) {
                    auto key_end_pos = line.find(':');
                    auto value_start_pos = line.find_first_not_of('\t', key_end_pos + 1);
                    auto key = line.substr(0, key_end_pos);
                    if (key == "Pid") {
                        stat.pid = extract_number(line, value_start_pos);
                    }
                    if (key == "VmPeak") {
                        stat.VmPeak = extract_number(line, value_start_pos, true);
                    }
                    if (key == "VmSize") {
                        stat.VmSize = extract_number(line, value_start_pos, true);
                    }
                    if (key == "VmHWM") {
                        stat.VmHWM = extract_number(line, value_start_pos, true);
                    }
                    if (key == "VmRSS") {
                        stat.VmRSS = extract_number(line, value_start_pos, true);
                    }
                    if (key == "VmData") {
                        stat.VmData = extract_number(line, value_start_pos, true);
                    }
                    if (key == "VmPTE") {
                        stat.VmPTE = extract_number(line, value_start_pos, true);
                    }
                }
            }

            return stat;
        }

};

