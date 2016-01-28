#pragma once

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <regex>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "timings.hpp"
#include "sdsl/int_vector.hpp"

typedef std::vector<uint64_t> string_type;

struct gapped_pattern {
    std::string raw_regexp;
    sdsl::int_vector<0> sdsl_regexp;
    std::vector<string_type> subpatterns;
    std::vector<std::string> gap_strs;
    std::vector<std::pair<uint64_t,uint64_t>> gaps;
    gapped_pattern(const std::string& p, bool string_patterns) : raw_regexp(p)
    {
        auto parse_subpattern = [&](std::string s) {
            if (string_patterns) {
                return string_type(s.begin(), s.end());
            } else {
                std::istringstream symbol_stream(s);
                string_type subpattern;
                uint64_t symbol;
                while (symbol_stream >> symbol)
                    subpattern.push_back(symbol);
                return subpattern;
            }
        };

        sdsl_regexp.resize(raw_regexp.size());
        for (size_t i=0; i<raw_regexp.size(); i++) {
            sdsl_regexp[i] = raw_regexp[i];
        }
        std::regex gap_regex(R"(\.\{[^,]*,[^,]*\})");
        auto gaps_begin = std::sregex_iterator(raw_regexp.begin(),raw_regexp.end(),gap_regex);
        auto gaps_end = std::sregex_iterator();
        int64_t last_gap_end = -1;
        for (auto i = gaps_begin; i != gaps_end; ++i) {
            std::smatch match = *i;
            std::string gap_str = match.str();
            gap_strs.push_back(gap_str);
            /* try parsing the gap  description */
            auto first_num_sep = gap_str.find(",");
            auto first_num_str = gap_str.substr(2,first_num_sep-2);
            auto second_num_str = gap_str.substr(first_num_sep+1);
            second_num_str.pop_back();
            uint64_t first_gap_num = std::stoull(first_num_str);
            uint64_t second_gap_num = std::stoull(second_num_str);
            if (first_gap_num > second_gap_num) {
                throw std::runtime_error("invalid gap description");
            }
            gaps.emplace_back(first_gap_num,second_gap_num);
            auto gap_pos = match.position();
            auto subptrlen = gap_pos - (last_gap_end+1);
            auto subpattern = raw_regexp.substr(last_gap_end+1,subptrlen);
            subpatterns.push_back(parse_subpattern(subpattern));
            last_gap_end = gap_pos + gap_str.size() - 1;
        }
        auto last_subpattern = raw_regexp.substr(last_gap_end+1);
        subpatterns.push_back(parse_subpattern(last_subpattern));
        LOG(INFO) << "PARSED('" << raw_regexp << "') -> GAPS = " << gaps << " SUBPATTERNS = " << subpatterns;
    };
};

struct gapped_search_result {
    std::vector<uint64_t> positions;
    gapped_search_result() = default;
    gapped_search_result(size_t n)
    {
        positions.resize(n);
    }
};

namespace utils
{

std::vector<gapped_pattern>
parse_pattern_file(std::string file_name, bool string_patterns)
{
    gm_timer tm("READ PATTERNS",true);
    std::vector<gapped_pattern> patterns;
    std::ifstream in(file_name);
    if (in) {
        std::string line;
        while (std::getline(in,line)) {
            try {
                gapped_pattern pat(line, string_patterns);
                patterns.push_back(pat);
            } catch (...) {
                LOG(ERROR) << "Could not parse pattern '" << line << "'. Skipped";
            }
        }
    } else {
        LOG(FATAL) << "Cannot open pattern file '" << file_name << "'";
    }
    LOG(INFO) << "read " << patterns.size() << " patterns";
    return patterns;
}

bool directory_exists(std::string dir)
{
    struct stat sb;
    const char* pathname = dir.c_str();
    if (stat(pathname, &sb) == 0 && (S_IFDIR & sb.st_mode)) {
        return true;
    }
    return false;
}

bool file_exists(std::string file_name)
{
    std::ifstream in(file_name);
    if (in) {
        in.close();
        return true;
    }
    return false;
}

void create_directory(std::string dir)
{
    if (!directory_exists(dir)) {
        if (mkdir(dir.c_str(), 0777) == -1) {
            LOG(FATAL) << "could not create directory";
        }
    }
}


}
