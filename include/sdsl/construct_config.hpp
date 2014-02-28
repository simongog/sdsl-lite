#ifndef INCLUDED_SDSL_CONSTRUCT_CONFIG
#define INCLUDED_SDSL_CONSTRUCT_CONFIG

#include "config.hpp"

namespace sdsl
{

class construct_config
{
    public:
        static byte_sa_algo_type byte_algo_sa;

        construct_config() = delete;
};

}

#endif
