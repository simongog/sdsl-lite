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

#include <sys/types.h> // for file_size
#include <sys/stat.h>  // for file_size
#include <iomanip>
#include <vector>
#include <string>

#include <type_traits>
#include <typeinfo>
#ifndef MSVC_COMPILER
#include <cxxabi.h>
#endif

namespace sdsl
{

namespace util
{

uint64_t _id_helper::id = 0;

std::string basename(std::string file)
{
    file = disk_file_name(file); // remove RAM-prefix
#ifdef MSVC_COMPILER
    char* c = _strdup((const char*)file.c_str());
    char file_name[_MAX_FNAME] = { 0 };
    ::_splitpath_s(c, NULL, 0, NULL, NULL, file_name, _MAX_FNAME, NULL, 0);
    std::string res(file_name);
#else
    char* c = strdup((const char*)file.c_str());
    std::string res = std::string(::basename(c));
#endif
    free(c);
    return res;
}

std::string dirname(std::string file)
{
    bool ram_file = is_ram_file(file);
    file = disk_file_name(file); // remove RAM-prefix
#ifdef MSVC_COMPILER
    char* c = _strdup((const char*)file.c_str());
    char dir_name[_MAX_DIR] = { 0 };
    char drive[_MAX_DRIVE] = {0};
    ::_splitpath_s(c, drive, _MAX_DRIVE, dir_name, _MAX_DIR, NULL,0, NULL,0);
    std::string res = std::string(drive) + std::string(dir_name);
#else
    char* c = strdup((const char*)file.c_str());
    std::string res = std::string(::dirname(c));
#endif
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
#ifdef MSVC_COMPILER
    return _getpid();
#else
    return getpid();
#endif
}

char* str_from_errno()
{
#ifdef MSVC_COMPILER
#pragma warning(disable:4996)
    return strerror(errno);
#pragma warning(default:4996)
#else
    return strerror(errno);
#endif
}


uint64_t id()
{
    return _id_helper::getId();
}

std::string demangle(const std::string& name)
{
#ifdef HAVE_CXA_DEMANGLE
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
    for (auto file_pair : file_map) {
        sdsl::remove(file_pair.second);
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

size_t file_size(const std::string& file)
{
    if (is_ram_file(file)) {
        return ram_fs::file_size(file);
    } else {
        struct stat fs;
        stat(file.c_str(), &fs);
        return fs.st_size;
    }
}

}// end namespace util

}// end namespace sdsl

