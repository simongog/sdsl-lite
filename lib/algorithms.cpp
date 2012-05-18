/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
#include "sdsl/algorithms.hpp"
#include <cstring>
#include <cmath>
#include <map>
#include <string>

typedef std::map<std::string, std::string> tMSS;

namespace sdsl
{

double algorithm::H_0(const unsigned char* c)
{
    size_t n = strlen((const char*)c);
    if (n==0)
        return 0;
    size_t cnt[256] = {};
    double res = 0;
    for (size_t i=0; i<n; ++i)
        ++cnt[c[i]];
    for (uint16_t i=1; i<256; ++i)
        if (cnt[i]>0) {
            res += (double)cnt[i]/(double)n * log2((double)n/(double)cnt[i]);
//			std::cerr<<"cnt["<<i<<"] = "<<cnt[i]<<"  n = "<<n<<std::endl;
        }
    return res;
}

double algorithm::H_0s(const unsigned char* c)
{
    size_t n = strlen((const char*)c);
    if (n==0)
        return 0;
    size_t cnt[256] = {};
    double res = 0;
    for (size_t i=0; i<n; ++i)
        ++cnt[c[i]];
    for (uint16_t i=1; i<256; ++i)
        if (cnt[i]>0) {
            res += (double)cnt[i]/(double)n * log2((double)n/(double)cnt[i]);
//			std::cerr<<"cnt["<<i<<"] = "<<cnt[i]<<"  n = "<<n<<std::endl;
        }
    if (res==0) { //only one character occures in the string
        return (floor(log2((double)n))+1)/(double)n;
    }
    return res;

}

}
