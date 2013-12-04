#ifndef SDSL_INC_RASTER_IMG
#define SDSL_INC_RASTER_IMG

#include <sdsl/bit_vectors.hpp>

namespace sdsl
{

struct raster_img {
    typedef uint64_t size_type;

    uint64_t     max_x; // max x value
    uint64_t     max_y; // max y value
    uint64_t     max_z; // max z value in the compacted range
    uint32_t     offset;
    bit_vector   value_map;
    int_vector<> data;

    //! Serializes the data structure into the given ostream
    uint64_t serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

    //! Loads the data structure from the given istream.
    void load(std::istream& in);
};

}

#endif
