#pragma once

#pragma once

#include <stdexcept>
#include <chrono>
#include <sdsl/qsufsort.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_mapper.hpp>

#include "utils.hpp"
#include "logging.hpp"
#include "constants.hpp"
#include "timings.hpp"


std::vector<std::string> collection_keys = { consts::KEY_TEXT };

struct collection {
    std::string path;
    std::map<std::string, std::string> file_map;
    collection() = default;
    collection(const std::string& p)
        : path(p + "/")
    {
        if (!utils::directory_exists(path)) {
            LOG(FATAL) << "collection path not found.";
            throw std::runtime_error("collection path not found.");
        }
        /* make sure the necessary files are present */
        if (!utils::file_exists(path + "/" + consts::KEY_PREFIX + consts::KEY_TEXT)) {
            LOG(FATAL) << "collection path does not contain text.";
            throw std::runtime_error("collection path does not contain text.");
        }
        /* register files that are present */
        for (const auto& key : collection_keys) {
            auto file_path = path + "/" + consts::KEY_PREFIX + key;
            if (utils::file_exists(file_path)) {
                file_map[key] = file_path;
                LOG(INFO) << "FOUND '" << key << "' at '" << file_path << "'";
            }
        }
    }
};
