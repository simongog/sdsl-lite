/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

#include "sdsl/util.hpp"
#include "cxxabi.h"
#include <sys/types.h> // for file_size
#include <sys/stat.h>  // for file_size
#include <iomanip>
#include <vector>
#include <string>

namespace sdsl
{

namespace util
{

uint64_t _id_helper::id = 0;
timeval stop_watch::m_first_t = {0,0};
rusage stop_watch::m_first_r = {{0,0},{0,0}};



std::string basename(std::string file)
{
    file = disk_file_name(file); // remove RAM-prefix
    char* c = strdup((const char*)file.c_str());
    std::string res = std::string(::basename(c));
    free(c);
    return res;
}

std::string dirname(std::string file)
{
    bool ram_file = is_ram_file(file);
    file = disk_file_name(file); // remove RAM-prefix
    char* c = strdup((const char*)file.c_str());
    std::string res = std::string(::dirname(c));
    free(c);
    if (ram_file) {
        if ("." == res) {
            res = ram_file_name("");
        } else if ("/" ==res) {
            res = ram_file_name(res);
        }
    }
    return res;
}

uint64_t pid()
{
    return getpid();
}

uint64_t id()
{
    return _id_helper::getId();
}

std::string demangle(const std::string& name)
{
#ifndef HAVE_CXA_DEMANGLE
    char buf[4096];
    size_t size = 4096;
    int status = 0;
    abi::__cxa_demangle(name.c_str(), buf, &size, &status);
    if (status==0)
        return std::string(buf);
    return name;
#else
    return name;
#endif
}

std::string demangle2(const std::string& name)
{
    std::string result = demangle(name);
    std::vector<std::string> words_to_delete;
    words_to_delete.push_back("sdsl::");
    words_to_delete.push_back("(unsigned char)");
    words_to_delete.push_back(", unsigned long");

    for (size_t k=0; k<words_to_delete.size(); ++k) {
        std::string w = words_to_delete[k];
        for (size_t i = result.find(w); i != std::string::npos; i = result.find(w, i)) {
            result.erase(i, w.length());
            ++i;
        }
    }
    size_t index = 0;
    std::string to_replace = "int_vector<1>";
    while ((index = result.find(to_replace, index)) != std::string::npos) {
        result.replace(index, to_replace.size(), "bit_vector");
    }
    return result;
}

void delete_all_files(tMSS& file_map)
{
    for (tMSS::iterator file_it=file_map.begin(); file_it!=file_map.end(); ++file_it) {
        sdsl::remove(file_it->second);
    }
    file_map.clear();
}

std::string to_latex_string(unsigned char c)
{
    if (c == '_')
        return "\\_";
    else if (c == '\0')
        return "\\$";
    else
        return to_string(c);
}

void set_verbose()
{
    verbose = true;
}

off_t file_size(const std::string& file)
{
    if (is_ram_file(file)) {
        return ram_fs::file_size(file);
    } else {
        struct stat filestatus;
        stat(file.c_str(), &filestatus);
        return filestatus.st_size;
    }
}

stop_watch::stop_watch() : m_ruse1(), m_ruse2(), m_timeOfDay1(), m_timeOfDay2()
{
    timeval t;
    t.tv_sec = 0; t.tv_usec = 0;
    m_ruse1.ru_utime = t; m_ruse1.ru_stime = t; // init m_ruse1
    m_ruse2.ru_utime = t; m_ruse2.ru_stime = t; // init m_ruse2
    m_timeOfDay1 = t; m_timeOfDay2 = t;
    if (m_first_t.tv_sec == 0) {
        gettimeofday(&m_first_t, 0);
    }
    if (m_first_r.ru_utime.tv_sec == 0 and m_first_r.ru_utime.tv_usec ==0) {
        getrusage(RUSAGE_SELF, &m_first_r);
    }
}

void stop_watch::start()
{
    gettimeofday(&m_timeOfDay1, 0);
    getrusage(RUSAGE_SELF, &m_ruse1);
}

void stop_watch::stop()
{
    getrusage(RUSAGE_SELF, &m_ruse2);
    gettimeofday(&m_timeOfDay2, 0);
}

double stop_watch::user_time()
{
    timeval t1, t2;
    t1 = m_ruse1.ru_utime;
    t2 = m_ruse2.ru_utime;
    return ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec)))/1000.0;
}

double stop_watch::sys_time()
{
    timeval t1, t2;
    t1 = m_ruse1.ru_stime;
    t2 = m_ruse2.ru_stime;
    return ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec)))/1000.0;
}

double stop_watch::real_time()
{
    double result = ((double)((m_timeOfDay2.tv_sec*1000000 + m_timeOfDay2.tv_usec)-(m_timeOfDay1.tv_sec*1000000 + m_timeOfDay1.tv_usec)))/1000.0;
    if (result < sys_time() + user_time())
        return sys_time()+user_time();
    return result;
}

uint64_t stop_watch::abs_real_time()
{
    uint64_t result = (((m_timeOfDay2.tv_sec*1000000 + m_timeOfDay2.tv_usec - (m_first_t.tv_sec*1000000 + m_first_t.tv_usec))))/1000;
    return result;
}

uint64_t stop_watch::abs_user_time()
{
    timeval t1, t2;
    t1 = m_first_r.ru_utime;
    t2 = m_ruse2.ru_utime;
    return (t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec))/1000;
}


uint64_t stop_watch::abs_sys_time()
{
    timeval t1, t2;
    t1 = m_first_r.ru_stime;
    t2 = m_ruse2.ru_stime;
    return (t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec))/1000;
}

uint64_t stop_watch::abs_page_faults()
{
    return m_ruse2.ru_majflt - m_first_r.ru_majflt; // does not work on my platform
}

std::string time_string()
{
    time_t rawtime;
    struct tm* timeinfo;
    char buffer[1024];
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 1024, "%Y-%m-%d-%H%M%S", timeinfo);
    return buffer;
}

}// end namespace util

}// end namespace sdsl

