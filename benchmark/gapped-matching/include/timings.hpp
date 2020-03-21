#pragma once

#include <chrono>
#include <numeric>

#include "logging.hpp"

using namespace std::chrono;
using watch = std::chrono::high_resolution_clock;

struct gm_timer {
    watch::time_point start;
    std::string name;
    bool output;
    gm_timer(const std::string& _n,bool o = false) : name(_n), output(o)
    {
        if (output) LOG(INFO) << "START(" << name << ")";
        start = watch::now();
    }
    ~gm_timer()
    {
        auto stop = watch::now();
        auto time_spent = stop-start;
        if (output) LOG(INFO) << "STOP(" << name << ") - " << duration_cast<milliseconds>(time_spent).count() / 1000.0f << " sec";
    }
    watch::duration
    elapsed() const
    {
        return watch::now() - start;
    }
};

struct timings_summary {
    watch::duration total;
    watch::duration mean;
    watch::duration max;
    watch::duration min;
    watch::duration median;
    watch::duration qrt_1st;
    watch::duration qrt_3rd;
};

struct timing_results {
    std::vector<watch::duration> timings;
    void add_timing(watch::duration d)
    {
        timings.emplace_back(d);
    }
    timings_summary summary()
    {
        timings_summary t;
        std::sort(timings.begin(),timings.end());
        t.min = timings.front();
        t.max = timings.back();
        t.median = timings[timings.size()/2];
        t.total = std::accumulate(timings.begin(),timings.end(),watch::duration(0));
        t.mean = t.total / timings.size();
        t.qrt_1st = timings[timings.size()/4];
        t.qrt_3rd = timings[timings.size()*3/4];
        return t;
    }
};
