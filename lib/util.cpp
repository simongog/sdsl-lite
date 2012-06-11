/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
#include "sdsl/structure_tree.hpp"
#include "cxxabi.h"
#include <vector>

namespace sdsl
{

namespace util
{

uint64_t _id_helper::id = 0;

std::string basename(const std::string& file_name)
{
    char* c = strdup((const char*)file_name.c_str());
    std::string res = std::string(::basename(c));
    free(c);
    return res;
}

std::string dirname(const std::string& file_name)
{
    char* c = strdup((const char*)file_name.c_str());
    std::string res = std::string(::dirname(c));
    free(c);
    return res;
}

uint64_t get_pid()
{
    return getpid();
}

std::string demangle(const char* name)
{
#ifndef HAVE_CXA_DEMANGLE
    char buf[4096];
    size_t size = 4096;
    int status = 0;
    abi::__cxa_demangle(name, buf, &size, &status);
    if (status==0)
        return std::string(buf);
    return std::string(name);
#else
    return std::string(name);
#endif
}

std::string demangle2(const char* name)
{
    std::string result = demangle(name);
    std::vector<std::string> words_to_delete;
    words_to_delete.push_back("sdsl::");
    words_to_delete.push_back("(unsigned char)");

    for (size_t k=0; k<words_to_delete.size(); ++k) {
        std::string w = words_to_delete[k];
        for (size_t i = result.find(w); i != std::string::npos; i = result.find(w, i)) {
            result.erase(i, w.length());
            ++i;
        }
    }
    return result;
}

void delete_all_files(tMSS& file_map)
{
    for (tMSS::iterator file_it=file_map.begin(); file_it!=file_map.end(); ++file_it) {
        std::remove(file_it->second.c_str());
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

template<>
size_t write_member<std::string>(const std::string& t, std::ostream& out, structure_tree_node* v, std::string name)
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(t));
    size_t written_bytes = 0;
    written_bytes += write_member(t.size(), out, child, "length");
    out.write(t.c_str(), t.size());
    written_bytes += t.size();
    structure_tree::add_size(v, written_bytes);
    return written_bytes;
}

template<>
void read_member<std::string>(std::string& t, std::istream& in)
{
    std::string::size_type size;
    read_member(size, in);
    char* buf = new char[size];
    in.read(buf, size);
    std::string temp(buf, size);
    delete [] buf;
    t.swap(temp);
}

template<>
bool load_from_file(void*& v, const char* file_name)
{
    return true;
}

bool load_from_file(char*& v, const char* file_name)
{
    if (v != NULL) {
        delete [] v;
        v = NULL;
    }
    std::ifstream in;
    in.open(file_name, std::ios::binary | std::ios::in);
    if (in) {
        const uint64_t SDSL_BLOCK_SIZE = (1<<20);
        uint64_t n=0, read = 0;
        char buf[SDSL_BLOCK_SIZE], *cp;
        do {
            in.read(buf, SDSL_BLOCK_SIZE);
            read = in.gcount();
            n+=read;
        } while (SDSL_BLOCK_SIZE == read);
        if (n==0)
            return false;
        v = new char[n+1];
        in.close();
        in.open(file_name);
        if (!in) {
            delete [] v;
            v = NULL;
            return false;
        }
        cp=v;
        do {
            in.read(cp, SDSL_BLOCK_SIZE);
            read = in.gcount();
            cp+= read;
        } while (SDSL_BLOCK_SIZE == read);
        *(v+n) = '\0';
        return true;
    } else
        return false;
}


}// end namespace util
}// end namespace sdsl

