#ifndef __FBB_WAVELET_TREE_H_INCLUDED
#define __FBB_WAVELET_TREE_H_INCLUDED

#include <cstdint>
#include <vector>
#include <queue>
#include <algorithm>

#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/io.hpp>


#define ADD_NAVIGATIONAL_BLOCK_HEADER  // very slightly increases space, but removes the navigational rank queries, default: ON
//#define SPARSE_SUPERBLOCK_MAPPING    // reduces the space at the cost of rank query slowdown, default: OFF
#define ALLOW_VARIABLE_BLOCK_SIZE      // permits using different block size in each superblock, reduces both time and space, default: ON


namespace utils
{

void compute_symbol_freq(const std::uint8_t* text, std::uint64_t text_length,
                         std::vector<std::uint64_t>& freq)
{
    std::fill(freq.begin(), freq.end(), 0UL);
    for (std::uint64_t i = 0; i < text_length; ++i)
        ++freq[text[i]];
}

void compute_huffman_code_lengths(std::vector<std::uint64_t>& freq,
                                  std::vector<std::uint64_t>& codelen)
{
    std::fill(codelen.begin(), codelen.end(), 0UL);
    typedef std::pair<std::uint64_t, std::vector<std::uint8_t> > pq_item;
    std::priority_queue<pq_item, std::vector<pq_item>, std::greater<pq_item> > pq;
    for (std::uint64_t i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
            std::vector<std::uint8_t> v { (std::uint8_t)i };
            pq.push(std::make_pair(freq[i], v));
        }
    }

    while (pq.size() > 1) {
        pq_item x = pq.top(); pq.pop();
        pq_item y = pq.top(); pq.pop();
        std::vector<std::uint8_t> v = x.second;
        v.insert(v.end(), y.second.begin(), y.second.end());
        for (std::uint8_t i : v) ++codelen[i];
        pq.push(std::make_pair(x.first + y.first, v));
    }
}

void assign_canonical_huffman_codes(std::vector<std::uint64_t>& freq,
                                    std::vector<std::uint64_t>& codelen, std::vector<std::uint64_t>& code)
{
    std::fill(code.begin(), code.end(), 0UL);
    std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0) sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
    std::sort(sym.begin(), sym.end());
    for (std::uint64_t c = 0, i = 0; i < sym.size(); ++i) {
        if (i == 0) code[sym[i].second] = c;
        else {
            c = (c + 1) << (sym[i].first - sym[i - 1].first);
            code[sym[i].second] = c;
        }
    }
}

void pre_encode_block(const std::uint8_t* block_ptr, std::uint64_t block_size,
                      std::uint64_t& block_sigma, std::uint64_t& block_tree_height,
                      std::uint64_t& block_bv_size, std::vector<uint16_t>& global_to_block_mapping)
{
    // Compute Huffman code lengths.
    std::vector<std::uint64_t> freq(256), codelen(256);
    compute_symbol_freq(block_ptr, block_size, freq);
    compute_huffman_code_lengths(freq, codelen);

    // Sort symbols by frequency.
    std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0) sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
    std::sort(sym.begin(), sym.end());

    // Fill in all fields.
    block_bv_size = 0;
    block_sigma = sym.size();
    block_tree_height = *std::max_element(codelen.begin(), codelen.end());
    for (std::uint64_t i = 0; i < sym.size(); ++i) {
        global_to_block_mapping[sym[i].second] = (std::uint16_t)i;
        if (sym.size() > 1)
            block_bv_size += freq[sym[i].second] * codelen[sym[i].second];
    }
}

void encode_block(const std::uint8_t* block_ptr, std::vector<std::uint64_t>& block_rank,
                  std::uint64_t block_size, sdsl::bit_vector& superblock_bv, std::uint64_t& ones_count,
                  std::uint64_t superblock_bv_offset, std::uint8_t* block_header_ptr)
{
    // Compute Huffman code.
    std::vector<std::uint64_t> freq(256), codelen(256), code(256);
    compute_symbol_freq(block_ptr, block_size, freq);
    compute_huffman_code_lengths(freq, codelen);
    assign_canonical_huffman_codes(freq, codelen, code);

    //----------------------------------------------------------------------------
    // NOTE: we identify the node by the number constructed as follows: the
    //   bits on the path from the root to the node are treated as bits of
    //   number (the first bit on the path is the most significant in the
    //   number). We then prepend this number (in binary) with 1, e.g., the
    //   node with path 011 is identified by number 10 = 1011. This way every
    //   possible node in the tree has a unique id. Moreover, sorted IDs
    //   correspond to nodes visited level-by-level (and left-to-right within
    //   level), which is the order in which we concatenate the bitvectors.
    //   It has other nice properties, e.g., the ID of sibling of node with
    //   ID x is (x^1) -- this comes in handy during bitvector construction.
    //----------------------------------------------------------------------------

    // Collect IDs of all internal nodes in the tree.
    std::uint64_t max_code_length = *std::max_element(codelen.begin(), codelen.end());
    std::vector<std::uint64_t> internal_node_ids;
    for (std::uint64_t i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
                std::uint64_t id = (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
                internal_node_ids.push_back(id);
            }
        }
    }

    // Compute the mapping from internal nodes to bitvectors (which are
    // numbered with consecuting numbers starting from 0, according to the
    // order in which they are concatenated).
    std::sort(internal_node_ids.begin(), internal_node_ids.end());
    internal_node_ids.erase(std::unique(internal_node_ids.begin(),
                                        internal_node_ids.end()), internal_node_ids.end());
    std::vector<std::uint64_t> internal_node_bv_id(1UL << max_code_length);
    for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
        internal_node_bv_id[internal_node_ids[i]] = i;

    // Compute the size of bitvector for every internal node.
    std::vector<std::uint64_t> internal_node_bv_size(internal_node_ids.size(), 0UL);
    for (std::uint64_t i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
                std::uint64_t id = (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
                internal_node_bv_size[internal_node_bv_id[id]] += freq[i];
            }
        }
    }

    // Sort symbols by frequency.
    std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0) sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
    std::sort(sym.begin(), sym.end());

    std::vector<std::uint64_t> ones_in_bv(internal_node_ids.size(), 0UL);  // used later
    ones_count = 0;
    if (sym.size() > 1) {
        // Allocate bitvectors for all internal nodes.
        std::vector<std::uint8_t>** internal_node_bv = new std::vector<std::uint8_t>* [internal_node_ids.size()];
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
            internal_node_bv[i] = new std::vector<std::uint8_t>(internal_node_bv_size[i], 0);

        // Fill in the bitvectors for all internal nodes.
        std::vector<std::uint64_t> node_visit_count(1UL << (max_code_length + 1), 0UL);
        for (std::uint64_t i = 0; i < block_size; ++i) {
            std::uint8_t sym = block_ptr[i];
            std::uint64_t pos = i;
            for (std::uint64_t depth = 0; depth < codelen[sym]; ++depth) {
                std::uint64_t id = (((1UL << codelen[sym]) | code[sym]) >> (codelen[sym] - depth));
                if (depth > 0) {
                    pos -= node_visit_count[id ^ 1];
                    ++node_visit_count[id];
                }
                if (code[sym] & (1UL << (codelen[sym] - depth - 1))) {
                    (*internal_node_bv[internal_node_bv_id[id]])[pos] = 1;
                    ones_in_bv[internal_node_bv_id[id]] += 1;
                    ++ones_count;
                }
            }
            ++node_visit_count[(1UL << codelen[sym]) | code[sym]];
        }

        // Append bitvectors of internal nodes to superblock bitvector.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
            for (std::uint64_t j = 0; j < internal_node_bv[i]->size(); ++j)
                superblock_bv[superblock_bv_offset++] = (*internal_node_bv[i])[j];

        // Clean up.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
            delete internal_node_bv[i];
        delete[] internal_node_bv;
    }

    //----------------------------------------------------------------------------
    // NOTE: level_total_freq[d] is the total frequency of symbols that have
    //   code length longer than d. In other words it's the total length of
    //   bitvectors associated with internal nodes at depth d in the tree.
    //   The fact, that above we use 2 bytes to ancode level_total_freq
    //   limits the block size to 2^16. The fact that we used 3 bytes to
    //   store rank at block boundary limits the superblocks to 2^24.
    //----------------------------------------------------------------------------

    // Fill in (the variable-size component of) the block header.
    std::vector<std::uint64_t> codelen_freq(max_code_length, 0UL);
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0 && codelen[i] < max_code_length)
            ++codelen_freq[codelen[i]];  // encoded with 1 byte
    std::vector<std::uint64_t> level_total_freq(max_code_length, 0UL);
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0)
            for (std::uint64_t depth = 1; depth < codelen[i]; ++depth)
                level_total_freq[depth] += freq[i];  // encoded with 2 bytes
    // Store the number of leaves and total size of bitvectors corresponding to internal
    // nodes at each level in the tree (except root level and deepest levels) minus one.
    for (std::uint64_t depth = 1; depth < max_code_length; ++depth) {
        *(block_header_ptr++) = (std::uint8_t)codelen_freq[depth];                      // leaf count, 1 byte
        *((std::uint16_t*)block_header_ptr) = (std::uint16_t)level_total_freq[depth] - 1;   // total bv size minus one, 2 bytes
        block_header_ptr += 2;
    }
    // Store rank value (wrt to superblock boundary) and global symbol for each leaf.
    for (std::uint64_t i = 0; i < sym.size(); ++i) {
        std::uint64_t symbol = (std::uint64_t)sym[i].second;  // encoded with 1 byte
        std::uint64_t rank_value = block_rank[symbol];        // encoded with 3 bytes
        *((std::uint32_t*)block_header_ptr) = (std::uint32_t)(symbol | (rank_value << 8));
        block_header_ptr += 4;
    }
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
    // For every internal node, store the number of 1-bits in the bitvector
    // corresponding to that node and all its left siblings (excl. leaves).
    std::uint64_t number_of_internal_nodes_current_level = 1;
    for (std::uint64_t depth = 0, ptr = 0; depth < max_code_length; ++depth) {
        std::uint64_t one_bits_current_level_count = 0;
        for (std::uint64_t j = 0; j < number_of_internal_nodes_current_level; ++j) {
            one_bits_current_level_count += ones_in_bv[ptr++];
            *((std::uint16_t*)block_header_ptr) = (std::uint16_t)one_bits_current_level_count;   // 2 bytes
            block_header_ptr += 2;
        }
        if (depth + 1 != max_code_length) {
            std::uint64_t next_level_leaf_count = codelen_freq[depth + 1];
            number_of_internal_nodes_current_level <<= 1;
            number_of_internal_nodes_current_level -= next_level_leaf_count;
        }
    }
#endif
}

void encode_block2(const std::uint8_t* block_ptr, std::uint64_t block_size,
                   sdsl::bit_vector& superblock_bv, std::uint64_t superblock_bv_offset)
{
    // Compute Huffman code.
    std::vector<std::uint64_t> freq(256), codelen(256), code(256);
    compute_symbol_freq(block_ptr, block_size, freq);
    compute_huffman_code_lengths(freq, codelen);
    assign_canonical_huffman_codes(freq, codelen, code);

    // Collect IDs of all internal nodes in the tree.
    std::uint64_t max_code_length = *std::max_element(codelen.begin(), codelen.end());
    std::vector<std::uint64_t> internal_node_ids;
    for (std::uint64_t i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
                std::uint64_t id = (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
                internal_node_ids.push_back(id);
            }
        }
    }

    // Compute the mapping from internal nodes to bitvectors (which are
    // numbered with consecuting numbers starting from 0, according to the
    // order in which they are concatenated).
    std::sort(internal_node_ids.begin(), internal_node_ids.end());
    internal_node_ids.erase(std::unique(internal_node_ids.begin(),
                                        internal_node_ids.end()), internal_node_ids.end());
    std::vector<std::uint64_t> internal_node_bv_id(1UL << max_code_length);
    for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
        internal_node_bv_id[internal_node_ids[i]] = i;

    // Compute the size of bitvector for every internal node.
    std::vector<std::uint64_t> internal_node_bv_size(internal_node_ids.size(), 0UL);
    for (std::uint64_t i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
            for (std::uint64_t depth = 0; depth < codelen[i]; ++depth) {
                std::uint64_t id = (((1UL << codelen[i]) | code[i]) >> (codelen[i] - depth));
                internal_node_bv_size[internal_node_bv_id[id]] += freq[i];
            }
        }
    }

    // Sort symbols by frequency.
    std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0) sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
    std::sort(sym.begin(), sym.end());

    std::vector<std::uint64_t> ones_in_bv(internal_node_ids.size(), 0UL);  // used later
    if (sym.size() > 1) {
        // Allocate bitvectors for all internal nodes.
        std::vector<std::uint8_t>** internal_node_bv = new std::vector<std::uint8_t>* [internal_node_ids.size()];
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
            internal_node_bv[i] = new std::vector<std::uint8_t>(internal_node_bv_size[i], 0);

        // Fill in the bitvectors for all internal nodes.
        std::vector<std::uint64_t> node_visit_count(1UL << (max_code_length + 1), 0UL);
        for (std::uint64_t i = 0; i < block_size; ++i) {
            std::uint8_t sym = block_ptr[i];
            std::uint64_t pos = i;
            for (std::uint64_t depth = 0; depth < codelen[sym]; ++depth) {
                std::uint64_t id = (((1UL << codelen[sym]) | code[sym]) >> (codelen[sym] - depth));
                if (depth > 0) {
                    pos -= node_visit_count[id ^ 1];
                    ++node_visit_count[id];
                }
                if (code[sym] & (1UL << (codelen[sym] - depth - 1))) {
                    (*internal_node_bv[internal_node_bv_id[id]])[pos] = 1;
                    ones_in_bv[internal_node_bv_id[id]] += 1;
                }
            }
            ++node_visit_count[(1UL << codelen[sym]) | code[sym]];
        }

        // Append bitvectors of internal nodes to superblock bitvector.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
            for (std::uint64_t j = 0; j < internal_node_bv[i]->size(); ++j)
                superblock_bv[superblock_bv_offset++] = (*internal_node_bv[i])[j];

        // Clean up.
        for (std::uint64_t i = 0; i < internal_node_ids.size(); ++i)
            delete internal_node_bv[i];
        delete[] internal_node_bv;
    }
}


//==============================================================================
// Restore the Huffman code of the given symbol in block alphabet (`block_c')
// from the block header and return via variables `code' and `codelen'. We
// assume tree_height > 0.
//==============================================================================
void restore_code_from_block_header(std::uint64_t block_c,
                                    const std::uint8_t* block_header_ptr, std::uint64_t tree_height,
                                    std::uint64_t& code, std::uint64_t& codelen)
{
    code = 0;
    codelen = 1;
    std::uint64_t leaf_count = 0;
    while (codelen < tree_height) {
        code <<= 1;
        std::uint64_t this_level_leaf_count = *block_header_ptr;
        if (leaf_count + this_level_leaf_count > block_c) {
            code += (block_c - leaf_count);
            break;
        } else {
            code += this_level_leaf_count;
            ++codelen;
            leaf_count += this_level_leaf_count;
            block_header_ptr += 3;
        }
    }
    if (codelen == tree_height) {
        code <<= 1;
        code += (block_c - leaf_count);
    }
}

// This function is for debugging purpose.
long double avg_wtree_traversal_len(const std::uint8_t* text, std::uint64_t text_length)
{
    // Compute Huffman code lengths.
    std::vector<std::uint64_t> freq(256), codelen(256);
    utils::compute_symbol_freq(text, text_length, freq);
    utils::compute_huffman_code_lengths(freq, codelen);

    // Sort symbols by frequency.
    std::vector<std::pair<std::uint64_t, std::uint8_t> > sym;
    for (std::uint64_t i = 0; i < 256; ++i)
        if (freq[i] > 0) sym.push_back(std::make_pair(codelen[i], (std::uint8_t)i));
    std::sort(sym.begin(), sym.end());

    // Compute result.
    long double result = 0;
    for (std::uint64_t i = 0; i < sym.size(); ++i)
        result += freq[sym[i].second] * codelen[sym[i].second];
    result /= (long double)text_length;

    return result;
}

}  // namespace utils


template<class t_bitvector       = sdsl::hyb_vector<>,
         class t_rank            = typename t_bitvector::rank_1_type,
#ifndef ALLOW_VARIABLE_BLOCK_SIZE
         std::uint64_t t_bs_log  = 14,
#endif
         std::uint64_t t_sbs_log = 20
         >
class wt_fbb
{
    public:
        typedef sdsl::wt_tag index_category;
        typedef std::uint64_t size_type;
        using alphabet_category = sdsl::byte_alphabet_tag;

    private:
#ifndef ALLOW_VARIABLE_BLOCK_SIZE
        static const std::uint64_t k_block_size_log;
        static const std::uint64_t k_block_size;
#endif
        static const std::uint64_t k_superblock_size_log;
        static const std::uint64_t k_superblock_size;
        static const std::uint64_t k_hyperblock_size;

        struct block_header_item {
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
            std::uint32_t m_bv_rank;                 // rank at the beginning of block vector
#endif
            std::uint32_t m_bv_offset;               // offset in the superblock bitvector
            std::uint32_t m_var_size_header_offset;  // offset in the array storing variable-size components of block header
            std::uint8_t  m_sigma;                   // block alphabet size - 1
            std::uint8_t  m_tree_height;             // height of the wavelet tree for block
        } __attribute__((packed));

        struct superblock_header_item {
            std::uint8_t m_sigma;                           // superblock alphabet size - 1
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            std::uint8_t m_block_size_log;
#endif
            t_bitvector m_bitvector;
            t_rank m_rank_support;
            std::vector<std::uint8_t> m_var_block_headers;   // variable-size components of block header
            std::vector<block_header_item> m_block_headers;  // fixed-size components of block header
#ifdef SPARSE_SUPERBLOCK_MAPPING
            std::vector<std::uint8_t> m_mapping_data;
            sdsl::bit_vector_il<> m_mapping_bv;
            typename sdsl::bit_vector_il<>::rank_1_type m_mapping_bv_rank_support;
            typename sdsl::bit_vector_il<>::select_1_type m_mapping_bv_select_support;
            std::uint32_t m_one_bit_count;
#else
            std::vector<std::uint8_t> m_mapping;            // mapping from superblock alphabet to block alphabet
#endif
        };

        friend std::uint64_t serialize(const superblock_header_item& shi, std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "")
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(shi));
            std::uint64_t written_bytes = 0;
            written_bytes += sdsl::serialize(shi.m_sigma, out, child, "sigma");
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            written_bytes += sdsl::serialize(shi.m_block_size_log, out, child, "block_size_log");
#endif
            written_bytes += sdsl::serialize(shi.m_bitvector, out, child, "bitvector");
            written_bytes += sdsl::serialize(shi.m_rank_support, out, child, "rank_support");
            written_bytes += sdsl::serialize(shi.m_var_block_headers, out, child, "var_block_headers");
            written_bytes += sdsl::serialize(shi.m_block_headers, out, child, "block_headers");
#ifdef SPARSE_SUPERBLOCK_MAPPING
            written_bytes += sdsl::serialize(shi.m_mapping_data, out, child, "mapping_data");
            written_bytes += sdsl::serialize(shi.m_mapping_bv, out, child, "mapping_bv");
            written_bytes += sdsl::serialize(shi.m_mapping_bv_rank_support, out, child, "rank_support");
            written_bytes += sdsl::serialize(shi.m_mapping_bv_select_support, out, child, "select_support");
            written_bytes += sdsl::serialize(shi.m_one_bit_count, out, child, "one_bit_count");
#else
            written_bytes += sdsl::serialize(shi.m_mapping, out, child, "mapping");
#endif
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        friend void load(superblock_header_item& shi, std::istream& in)
        {
            sdsl::load(shi.m_sigma, in);
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            sdsl::load(shi.m_block_size_log, in);
#endif
            sdsl::load(shi.m_bitvector, in);
            sdsl::load(shi.m_rank_support, in);
            sdsl::load(shi.m_var_block_headers, in);
            sdsl::load(shi.m_block_headers, in);
#ifdef SPARSE_SUPERBLOCK_MAPPING
            sdsl::load(shi.m_mapping_data);
            sdsl::load(shi.m_mapping_bv);
            sdsl::load(shi.m_rank_support);
            sdsl::load(shi.m_select_support);
            sdsl::load(shi.m_one_bit_count);
            shi.m_rank_support.set_vector(&shi.m_mapping_bv);
            shi.m_select_support.set_vector(&shi.m_mapping_bv);
#else
            sdsl::load(shi.m_mapping, in);
#endif
        }

    public:
        void collect_sizes(std::uint64_t& bitvectors, std::uint64_t& fixed_block_headers, std::uint64_t& variable_block_headers,
                           std::uint64_t& superblock_mapping, std::uint64_t& superblock_headers, std::uint64_t& other)
        {
            std::uint64_t n_superblocks = (m_size + k_superblock_size - 1) / k_superblock_size;

            bitvectors = 0;
            for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id)
                bitvectors += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_bitvector) +
                              sdsl::size_in_bytes(m_superblock_header[superblock_id].m_rank_support);

            fixed_block_headers = 0;
            for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id)
                fixed_block_headers += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_block_headers);

            variable_block_headers = 0;
            for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id)
                variable_block_headers += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_var_block_headers);

            superblock_mapping = 0;
            for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id) {
#ifdef SPARSE_SUPERBLOCK_MAPPING
                superblock_mapping += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_mapping_data);
                superblock_mapping += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_mapping_bv);
                superblock_mapping += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_mapping_bv_rank_support);
                superblock_mapping += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_mapping_bv_select_support);
                superblock_mapping += sizeof(superblock_header_item::m_one_bit_count);
#else
                superblock_mapping += sdsl::size_in_bytes(m_superblock_header[superblock_id].m_mapping);
#endif
            }

            superblock_headers = 0;
            superblock_headers += sdsl::size_in_bytes(m_global_mapping);
            superblock_headers += sdsl::size_in_bytes(m_superblock_rank);
            superblock_headers += n_superblocks * sizeof(superblock_header_item::m_sigma);
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            superblock_headers += n_superblocks * sizeof(superblock_header_item::m_block_size_log);
#endif

            other += sdsl::size_in_bytes(m_count);
            other += sdsl::size_in_bytes(m_hyperblock_rank);
            other += sizeof(m_superblock_header.size());
            other += sizeof(m_size);
        }


    private:
        std::uint64_t m_size;                                             // original sequence length
        std::vector<std::uint64_t> m_count;                               // global symbol counts
        std::vector<std::uint64_t> m_hyperblock_rank;                     // ranks at hyperblock boundary
        std::vector<std::uint32_t> m_superblock_rank;                     // ranks at superblock boundary
        std::vector<std::uint8_t> m_global_mapping;                       // mapping from global alphabet to superblock alphabet
        std::vector<superblock_header_item> m_superblock_header;          // superblock headers

        void copy(const wt_fbb& tree)
        {
            m_size = tree.m_size;
            m_count = tree.m_count;
            m_hyperblock_rank = tree.m_hyperblock_rank;
            m_superblock_rank = tree.m_superblock_rank;
            m_global_mapping = tree.m_global_mapping;
            m_superblock_header = tree.m_superblock_header;
        }

    public:
        // Default constructor
        wt_fbb() = default;

        // Copy constructor
        wt_fbb(const wt_fbb& tree)
        {
            copy(tree);
        }

        // Move constructor
        wt_fbb(wt_fbb&& tree)
        {
            *this = std::move(tree);
        }

        // Constructor
        wt_fbb(sdsl::int_vector_buffer<(std::uint8_t)8>& text_buf, std::uint64_t text_length)
        {
            init(text_buf, text_length);
        }

        // Constructor
        wt_fbb(const std::uint8_t* text, std::uint64_t text_length)
        {
            init(text, text_length);
        }

    private:
        void encode_blocks_in_superblock(const std::uint8_t* superblock_ptr, std::uint64_t superblock_id, std::uint64_t block_size_log)
        {
            std::uint64_t block_size = (1UL << block_size_log);
            std::uint64_t superblock_beg = superblock_id * k_superblock_size;
            std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
            std::uint64_t this_superblock_size = superblock_end - superblock_beg;
            std::uint64_t this_superblock_sigma = (std::uint64_t)m_superblock_header[superblock_id].m_sigma + 1;
#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            m_superblock_header[superblock_id].m_block_size_log = block_size_log;
#endif

#ifdef SPARSE_SUPERBLOCK_MAPPING
            // Create remporary containers for superblock mapping.
            sdsl::bit_vector this_superblock_mapping_bv(this_superblock_sigma * (k_superblock_size / block_size), 0);
            std::vector<std::uint8_t> this_superblock_mapping_data[256];
            std::uint32_t one_bit_count = 0;
#else
            // Allocate the array storing mapping from superblock alphabet to block alpabets.
            m_superblock_header[superblock_id].m_mapping =
                std::vector<std::uint8_t>(this_superblock_sigma * (k_superblock_size / block_size), 255);
#endif

            // Fill in block/superblock headers and compute sizes of
            // variable-size component of block headers and bitvectors.
            std::uint64_t superblock_bv_size = 0;
            std::uint64_t variable_block_header_size = 0;
            std::uint64_t blocks_in_this_superblock = (this_superblock_size + block_size - 1) / block_size;
            m_superblock_header[superblock_id].m_block_headers = std::vector<block_header_item>(blocks_in_this_superblock);
            for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
                std::uint64_t block_beg = block_id * block_size;
                std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
                std::uint64_t this_block_size = block_end - block_beg;
                const std::uint8_t* block_ptr = superblock_ptr + block_beg;

                // Pre-encode the block.
                std::uint64_t tree_height = 0, bv_size = 0, sigma = 0;
                std::vector<std::uint16_t> global_to_block_mapping(256, 256);
                utils::pre_encode_block(block_ptr, this_block_size, sigma, tree_height, bv_size, global_to_block_mapping);

                // Fill in the block header.
                m_superblock_header[superblock_id].m_block_headers[block_id].m_bv_offset = (std::uint32_t)superblock_bv_size;
                m_superblock_header[superblock_id].m_block_headers[block_id].m_var_size_header_offset = (std::uint32_t)variable_block_header_size;
                m_superblock_header[superblock_id].m_block_headers[block_id].m_tree_height = (std::uint8_t)tree_height;
                m_superblock_header[superblock_id].m_block_headers[block_id].m_sigma = (std::uint8_t)(sigma - 1);

                // Compute sparse representation of superblock mapping.
                for (std::uint64_t i = 0; i < 256; ++i) {
                    if (global_to_block_mapping[i] != 256) {
                        std::uint8_t superblock_char = m_global_mapping[superblock_id * 256 + i];
                        std::uint64_t addr = (std::uint64_t)superblock_char * (k_superblock_size / block_size) + block_id;
                        std::uint8_t mapping_val = std::min((std::uint16_t)254, global_to_block_mapping[i]);
#ifdef SPARSE_SUPERBLOCK_MAPPING
                        this_superblock_mapping_bv[addr] = 1;
                        this_superblock_mapping_data[superblock_char].push_back(mapping_val);
                        ++one_bit_count;
#else
                        m_superblock_header[superblock_id].m_mapping[addr] = mapping_val;
#endif
                    }
                }

                // Update the size of superblock bitvector.
                superblock_bv_size += bv_size;

                // Update the total size of variable-size block headers.
                if (tree_height > 1)
                    variable_block_header_size += (tree_height - 1) * 3;  // info about each level in the tree
                variable_block_header_size += sigma * 4;                // info about each leaf in the tree
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
                variable_block_header_size += (sigma - 1) * 2;          // additional navigational info
#endif
            }

#ifdef SPARSE_SUPERBLOCK_MAPPING
            m_superblock_header[superblock_id].m_mapping_bv = sdsl::bit_vector_il<>(this_superblock_mapping_bv);
            m_superblock_header[superblock_id].m_mapping_bv_rank_support = sdsl::bit_vector_il<>::rank_1_type(&m_superblock_header[superblock_id].m_mapping_bv);
            m_superblock_header[superblock_id].m_mapping_bv_select_support = sdsl::bit_vector_il<>::select_1_type(&m_superblock_header[superblock_id].m_mapping_bv);

            std::uint64_t total_sparse_mapping_size = 0;
            for (std::uint64_t i = 0; i < 256; ++i)
                total_sparse_mapping_size += this_superblock_mapping_data[i].size();
            m_superblock_header[superblock_id].m_mapping_data = std::vector<std::uint8_t>(total_sparse_mapping_size);
            std::uint64_t mapping_data_ptr = 0;
            for (std::uint64_t i = 0; i < 256; ++i)
                for (std::uint64_t j = 0; j < this_superblock_mapping_data[i].size(); ++j)
                    m_superblock_header[superblock_id].m_mapping_data[mapping_data_ptr++] = this_superblock_mapping_data[i][j];
            m_superblock_header[superblock_id].m_one_bit_count = one_bit_count;
#endif

            // Allocate the variable size block header.
            m_superblock_header[superblock_id].m_var_block_headers = std::vector<std::uint8_t>(variable_block_header_size);

            // Fill in superblock bitvector and block headers.
            std::uint64_t bv_rank = 0;
            sdsl::bit_vector superblock_bv(superblock_bv_size);
            std::vector<std::uint64_t> block_rank(256, 0UL);
            for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
                std::uint64_t block_beg = block_id * block_size;
                std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
                std::uint64_t this_block_size = block_end - block_beg;
                const std::uint8_t* block_ptr = superblock_ptr + block_beg;

                // Fill in the variable size block header and block bitvector.
                std::uint64_t ones_count = 0;
                std::uint64_t superblock_bv_offset = m_superblock_header[superblock_id].m_block_headers[block_id].m_bv_offset;
                std::uint64_t variable_block_header_offset = m_superblock_header[superblock_id].m_block_headers[block_id].m_var_size_header_offset;
                utils::encode_block(block_ptr, block_rank, this_block_size, superblock_bv, ones_count,
                                    superblock_bv_offset, m_superblock_header[superblock_id].m_var_block_headers.data() + variable_block_header_offset);

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
                m_superblock_header[superblock_id].m_block_headers[block_id].m_bv_rank = bv_rank;
#endif

                // Update block rank.
                bv_rank += ones_count;
                for (std::uint64_t i = 0; i < this_block_size; ++i)
                    ++block_rank[block_ptr[i]];
            }

            // Convert the superblock bitvector to final encoding and store.
            m_superblock_header[superblock_id].m_bitvector = t_bitvector(superblock_bv);
            m_superblock_header[superblock_id].m_rank_support = t_rank(&m_superblock_header[superblock_id].m_bitvector);
        }

        std::uint64_t size_of_blocks_encoding(const std::uint8_t* superblock_ptr, std::uint64_t superblock_id, std::uint64_t block_size)
        {
            std::uint64_t result = 0;
            std::uint64_t superblock_beg = superblock_id * k_superblock_size;
            std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
            std::uint64_t this_superblock_size = superblock_end - superblock_beg;
            std::uint64_t this_superblock_sigma = (std::uint64_t)m_superblock_header[superblock_id].m_sigma + 1;

#ifdef SPARSE_SUPERBLOCK_MAPPING
            // Create remporary containers for superblock mapping.
            sdsl::bit_vector this_superblock_mapping_bv(this_superblock_sigma * (k_superblock_size / block_size), 0);
            std::vector<std::uint8_t> this_superblock_mapping_data[256];
#else
            // Allocate the array storing mapping from superblock alphabet to block alpabets.
            result += this_superblock_sigma * (k_superblock_size / block_size);
#endif

            // Fill in block/superblock headers and compute sizes of
            // variable-size component of block headers and bitvectors.
            std::uint64_t superblock_bv_size = 0;
            std::uint64_t variable_block_header_size = 0;
            std::uint64_t blocks_in_this_superblock = (this_superblock_size + block_size - 1) / block_size;
            result += sizeof(block_header_item) * blocks_in_this_superblock;

            std::vector<std::uint32_t> bv_offset(blocks_in_this_superblock);
            for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
                std::uint64_t block_beg = block_id * block_size;
                std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
                std::uint64_t this_block_size = block_end - block_beg;
                const std::uint8_t* block_ptr = superblock_ptr + block_beg;

                // Pre-encode the block.
                std::uint64_t tree_height = 0, bv_size = 0, sigma = 0;
                std::vector<std::uint16_t> global_to_block_mapping(256, 256);
                utils::pre_encode_block(block_ptr, this_block_size, sigma, tree_height, bv_size, global_to_block_mapping);

                // Fill in the block header.
                bv_offset[block_id] = (std::uint32_t)superblock_bv_size;

                // Compute sparse representation of superblock mapping.
#ifdef SPARSE_SUPERBLOCK_MAPPING
                for (std::uint64_t i = 0; i < 256; ++i) {
                    if (global_to_block_mapping[i] != 256) {
                        std::uint8_t superblock_char = m_global_mapping[superblock_id * 256 + i];
                        std::uint64_t addr = (std::uint64_t)superblock_char * (k_superblock_size / block_size) + block_id;
                        std::uint8_t mapping_val = std::min((std::uint16_t)254, global_to_block_mapping[i]);
                        this_superblock_mapping_bv[addr] = 1;
                        this_superblock_mapping_data[superblock_char].push_back(mapping_val);
                    }
                }
#endif

                // Update the size of superblock bitvector.
                superblock_bv_size += bv_size;

                // Update the total size of variable-size block headers.
                if (tree_height > 1)
                    variable_block_header_size += (tree_height - 1) * 3;  // info about each level in the tree
                variable_block_header_size += sigma * 4;                // info about each leaf in the tree
#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
                variable_block_header_size += (sigma - 1) * 2;          // additional navigational info
#endif
            }

#ifdef SPARSE_SUPERBLOCK_MAPPING
            sdsl::bit_vector_il<> temp_m_mapping_bv(this_superblock_mapping_bv);
            result += sdsl::size_in_bytes(temp_m_mapping_bv);
            for (std::uint64_t i = 0; i < 256; ++i)
                result += this_superblock_mapping_data[i].size();
#endif
            result += variable_block_header_size;

            sdsl::bit_vector superblock_bv(superblock_bv_size);
            for (std::uint64_t block_id = 0; block_id < blocks_in_this_superblock; ++block_id) {
                std::uint64_t block_beg = block_id * block_size;
                std::uint64_t block_end = std::min(block_beg + block_size, this_superblock_size);
                std::uint64_t this_block_size = block_end - block_beg;
                const std::uint8_t* block_ptr = superblock_ptr + block_beg;
                utils::encode_block2(block_ptr, this_block_size, superblock_bv, (std::uint64_t)bv_offset[block_id]);
            }
            t_bitvector temp_m_bitvector(superblock_bv);
            t_rank temp_m_rank_support(&temp_m_bitvector);
            result += sdsl::size_in_bytes(temp_m_bitvector);
            result += sdsl::size_in_bytes(temp_m_rank_support);

            return result;
        }

        void encode_superblock(const std::uint8_t* superblock_ptr, std::uint64_t superblock_id)
        {
            std::uint64_t superblock_beg = superblock_id * k_superblock_size;
            std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
            std::uint64_t this_superblock_size = superblock_end - superblock_beg;

            // Store ranks at hyperblock boundary.
            std::uint64_t hyperblock_id = (superblock_id * k_superblock_size) / k_hyperblock_size;
            if (!(superblock_id * k_superblock_size % k_hyperblock_size))
                for (std::uint64_t i = 0; i < 256; ++i)
                    m_hyperblock_rank[hyperblock_id * 256 + i] = m_count[i];

            // Store ranks at superblock boundary.
            for (std::uint64_t i = 0; i < 256; ++i)
                m_superblock_rank[superblock_id * 256 + i] = m_count[i] - m_hyperblock_rank[hyperblock_id * 256 + i];

            // Update symbol counts.
            for (std::uint64_t i = 0; i < this_superblock_size; ++i)
                ++m_count[(std::uint8_t)superblock_ptr[i]];

            // Compute superblock sigma and mapping from global alpabet to superblock alphabet.
            std::uint64_t this_superblock_sigma = 0;
            for (std::uint64_t i = 0; i < 256; ++i)
                if (m_superblock_rank[superblock_id * 256 + i] + m_hyperblock_rank[hyperblock_id * 256 + i] != m_count[i])
                    m_global_mapping[superblock_id * 256 + i] = this_superblock_sigma++;
            m_superblock_header[superblock_id].m_sigma = this_superblock_sigma - 1;

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            // Try few different block sizes.
            std::uint64_t best_block_size_log = 0;
            std::uint64_t best_encoding_size = (1UL << 60);
            for (std::uint64_t block_size_log = (std::uint64_t)std::max(0L, (std::int64_t)std::min(k_superblock_size_log, 16UL) - 7);
                 block_size_log <= std::min(k_superblock_size_log, 16UL); ++block_size_log) {
                std::uint64_t encoding_size = size_of_blocks_encoding(superblock_ptr, superblock_id, (1UL << block_size_log));
                if (encoding_size < best_encoding_size) {
                    best_block_size_log = block_size_log;
                    best_encoding_size = encoding_size;
                }
            }

            // Encode blocks inside superblock.
            encode_blocks_in_superblock(superblock_ptr, superblock_id, best_block_size_log);
#else
            encode_blocks_in_superblock(superblock_ptr, superblock_id, k_block_size_log);
#endif
        }

        void init(const std::uint8_t* text, std::uint64_t text_length)
        {
            m_size = text_length;
            m_count = std::vector<std::uint64_t>(256, 0UL);
            std::uint64_t n_superblocks = (m_size + k_superblock_size - 1) / k_superblock_size;
            std::uint64_t n_hyperblocks = (m_size + k_hyperblock_size - 1) / k_hyperblock_size;

            // Allocate headers.
            m_hyperblock_rank = std::vector<std::uint64_t>(n_hyperblocks * 256);
            m_superblock_rank = std::vector<std::uint32_t>(n_superblocks * 256);
            m_global_mapping = std::vector<std::uint8_t>(n_superblocks * 256, 255);
            m_superblock_header = std::vector<superblock_header_item>(n_superblocks);

            // Encode superblocks, left to right.
            for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id) {
                std::uint64_t superblock_beg = superblock_id * k_superblock_size;
                std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
                std::uint64_t this_superblock_size = superblock_end - superblock_beg;
                std::uint8_t* superblock_ptr = new std::uint8_t[this_superblock_size];
                std::copy(text + superblock_beg, text + superblock_end, superblock_ptr);

                encode_superblock(superblock_ptr, superblock_id);
                delete[] superblock_ptr;
            }
        }

        void init(sdsl::int_vector_buffer<(std::uint8_t)8>& text_buf, std::uint64_t text_length)
        {
            m_size = text_length;
            m_count = std::vector<std::uint64_t>(256, 0UL);
            std::uint64_t n_superblocks = (m_size + k_superblock_size - 1) / k_superblock_size;
            std::uint64_t n_hyperblocks = (m_size + k_hyperblock_size - 1) / k_hyperblock_size;

            // Allocate headers.
            m_hyperblock_rank = std::vector<std::uint64_t>(n_hyperblocks * 256);
            m_superblock_rank = std::vector<std::uint32_t>(n_superblocks * 256);
            m_global_mapping = std::vector<std::uint8_t>(n_superblocks * 256, 255);
            m_superblock_header = std::vector<superblock_header_item>(n_superblocks);

            // Encode superblocks, left to rightt.
            for (std::uint64_t superblock_id = 0; superblock_id < n_superblocks; ++superblock_id) {
                std::uint64_t superblock_beg = superblock_id * k_superblock_size;
                std::uint64_t superblock_end = std::min(superblock_beg + k_superblock_size, m_size);
                std::uint64_t this_superblock_size = superblock_end - superblock_beg;
                std::uint8_t* superblock_ptr = new std::uint8_t[this_superblock_size];
                for (std::uint64_t i = 0; i < this_superblock_size; ++i)
                    superblock_ptr[i] = text_buf[superblock_beg + i];

                encode_superblock(superblock_ptr, superblock_id);
                delete[] superblock_ptr;
            }
        }

    public:
        // We assume 0 <= i <= m_size.
        std::uint64_t rank(std::uint64_t i, std::uint8_t c) const
        {
            if (i == 0) return 0;
            else if (i == m_size) return m_count[c];

            std::uint64_t hyperblock_id = i / k_hyperblock_size;
            std::uint64_t superblock_id = i / k_superblock_size;
            std::uint64_t superblock_i = i % k_superblock_size;
            std::uint8_t superblock_c = m_global_mapping[superblock_id * 256 + c];
            std::uint64_t superblock_sigma = (std::uint64_t)m_superblock_header[superblock_id].m_sigma + 1;

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
            std::uint64_t block_size_log = m_superblock_header[superblock_id].m_block_size_log;
#else
            std::uint64_t block_size_log = k_block_size_log;
#endif

            std::uint64_t block_size = (1UL << block_size_log);
            std::uint64_t blocks_in_superblock_log = (k_superblock_size_log - block_size_log);
            std::uint64_t block_i = i & (block_size - 1);
            std::uint64_t this_block_size = std::min(block_size, m_size - (i - block_i));
            std::uint64_t block_id = (superblock_i >> block_size_log);
            std::uint64_t rank_at_superblock_boundary = m_superblock_rank[superblock_id * 256 + c];
            std::uint64_t rank_at_hyperblock_boundary = m_hyperblock_rank[hyperblock_id * 256 + c];
            if (superblock_c >= superblock_sigma)  // special case: c does not occur in the superblock
                return rank_at_hyperblock_boundary + rank_at_superblock_boundary;

#ifdef SPARSE_SUPERBLOCK_MAPPING
            std::uint64_t addr = ((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id;
            bool mapping_bit = m_superblock_header[superblock_id].m_mapping_bv[addr];
            std::uint64_t mapping_rank = m_superblock_header[superblock_id].m_mapping_bv_rank_support.rank(addr);
            std::uint8_t block_c = mapping_bit ? m_superblock_header[superblock_id].m_mapping_data[mapping_rank] : 255;
#else
            const std::uint8_t* superblock_mapping_ptr = m_superblock_header[superblock_id].m_mapping.data();
            std::uint8_t block_c = superblock_mapping_ptr[((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id];
#endif

            if (block_c == 255) {  // special case: c does not occur in the block
#ifdef SPARSE_SUPERBLOCK_MAPPING
                std::uint64_t select_ret = 0;
                if (mapping_rank != m_superblock_header[superblock_id].m_one_bit_count)
                    select_ret = m_superblock_header[superblock_id].m_mapping_bv_select_support.select(mapping_rank + 1);
                if (mapping_rank == m_superblock_header[superblock_id].m_one_bit_count ||
                    select_ret >= ((superblock_c + 1) << blocks_in_superblock_log)) {
#else
                // Find the closest block to the right in which c occurs.
                ++block_id;
                std::uint64_t blocks_in_superblock = (1UL << blocks_in_superblock_log);
                while (block_id < blocks_in_superblock &&
                       superblock_mapping_ptr[((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id] == 255)
                    ++block_id;
                if (block_id == blocks_in_superblock) {
#endif
                    // Return the answer from superblock header or global counts.
                    if ((superblock_id + 1) * k_superblock_size >= m_size) return m_count[c];
                    else return rank_at_hyperblock_boundary + m_superblock_rank[(superblock_id + 1) * 256 + c];
                } else {
#ifdef SPARSE_SUPERBLOCK_MAPPING
                    block_c = m_superblock_header[superblock_id].m_mapping_data[mapping_rank];
                    block_id = m_superblock_header[superblock_id].m_mapping_bv_select_support.select(mapping_rank + 1) -
                               ((std::uint64_t)superblock_c << blocks_in_superblock_log);
#else
                    block_c = superblock_mapping_ptr[((std::uint64_t)superblock_c << blocks_in_superblock_log) + block_id];
#endif

                    // Return the rank value from block header.
                    std::uint64_t var_size_header_offset = m_superblock_header[superblock_id].m_block_headers[block_id].m_var_size_header_offset;
                    std::uint64_t block_tree_height = m_superblock_header[superblock_id].m_block_headers[block_id].m_tree_height;
                    const std::uint8_t* variable_block_header_ptr = m_superblock_header[superblock_id].m_var_block_headers.data() + var_size_header_offset;
                    if (block_tree_height > 0)
                        variable_block_header_ptr += (block_tree_height - 1) * 3;
                    std::uint32_t* variable_block_header_ptr32 = (std::uint32_t*)variable_block_header_ptr;
                    if ((variable_block_header_ptr32[block_c] & 0xff) != c)  // special case: block_c == 254
                        ++block_c;
                    std::uint64_t rank_at_block_boundary = (variable_block_header_ptr32[block_c] >> 8);
                    return rank_at_hyperblock_boundary + rank_at_superblock_boundary + rank_at_block_boundary;
                }
            }

            // Compute rank at block boundary.
            std::uint64_t var_size_header_offset = m_superblock_header[superblock_id].m_block_headers[block_id].m_var_size_header_offset;
            std::uint64_t block_tree_height = m_superblock_header[superblock_id].m_block_headers[block_id].m_tree_height;
            const std::uint8_t* variable_block_header_ptr = m_superblock_header[superblock_id].m_var_block_headers.data() + var_size_header_offset;
            const std::uint8_t* variable_block_header_ptr_temp = variable_block_header_ptr;
            if (block_tree_height > 0)
                variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
            const std::uint32_t* variable_block_header_ptr32 = (std::uint32_t*)variable_block_header_ptr_temp;
            if ((variable_block_header_ptr32[block_c] & 255) != c)  // special case: block_c == 254
                ++block_c;
            std::uint64_t rank_at_block_boundary = (variable_block_header_ptr32[block_c] >> 8);

            // Answer rank query inside block.
            if (block_tree_height == 0) // special case, a block was a run of single symbol.
                return rank_at_hyperblock_boundary + rank_at_superblock_boundary + rank_at_block_boundary + block_i;
            std::uint64_t code = 0;
            std::uint64_t codelen = 0;
            utils::restore_code_from_block_header(block_c, variable_block_header_ptr, block_tree_height, code, codelen);

#ifdef ADD_NAVIGATIONAL_BLOCK_HEADER
            std::uint64_t bv_rank = m_superblock_header[superblock_id].m_block_headers[block_id].m_bv_rank;      // rank (in s-block bv) at the beginning of current level
            std::uint64_t bv_offset = m_superblock_header[superblock_id].m_block_headers[block_id].m_bv_offset;  // starting pos. (in s-block bv) of bv at current level
            std::uint64_t int_nodes_count = 1;                         // number of internal nodes at current level
            std::uint64_t left_siblings_count = 0;                     // number of left siblings (excl. leaves) of current node
            std::uint64_t left_siblings_total_bv_size = 0;             // total size of bv's corresponding to left siblings of current node (excl. leaves)
            std::uint64_t cur_node_bv_size = this_block_size;          // size of bitvector of current node
            std::uint64_t cur_depth_total_bv_size = cur_node_bv_size;  // total size of bitvectors at current depth
            std::uint64_t cur_node_rank = block_i;                     // the rank value that we refine at each level

            // We maintain second pointer to variable size block header. It is used to extract
            // total number of 1s in bitvectors corresponding to left siblings of current node
            // (excl. leaves) and also including 1s in the bitvector of current node.
            std::uint64_t block_sigma = (std::uint64_t)m_superblock_header[superblock_id].m_block_headers[block_id].m_sigma + 1;
            variable_block_header_ptr_temp = variable_block_header_ptr;
            if (block_tree_height > 0)
                variable_block_header_ptr_temp += (block_tree_height - 1) * 3;
            variable_block_header_ptr_temp += block_sigma * 4;
            std::uint16_t* second_variable_block_header_ptr = (std::uint16_t*)variable_block_header_ptr_temp;

            // Traverse the tree.
            for (std::uint64_t depth = 0; depth < codelen; ++depth) {
                // Compute the number of 1s in current node.
                std::uint64_t rank1 = m_superblock_header[superblock_id].m_rank_support.rank(bv_offset + left_siblings_total_bv_size + cur_node_rank);
                std::uint64_t left_siblings_total_ones_count = 0;
                if (left_siblings_count > 0)
                    left_siblings_total_ones_count = (std::uint64_t)second_variable_block_header_ptr[left_siblings_count - 1];
                rank1 -= bv_rank + left_siblings_total_ones_count;

                // Compute remaining stats about current node.
                std::uint64_t cur_node_one_count = (std::uint64_t)second_variable_block_header_ptr[left_siblings_count] - left_siblings_total_ones_count;
                std::uint64_t cur_node_zero_count = cur_node_bv_size - cur_node_one_count;
                std::uint64_t rank0 = cur_node_rank - rank1;

                // Update nagivational info.
                bv_rank += (std::uint64_t)second_variable_block_header_ptr[int_nodes_count - 1];
                second_variable_block_header_ptr += int_nodes_count;
                left_siblings_count <<= 1;

                // Update rank.
                std::uint64_t next_bit = (code & (1UL << (codelen - depth - 1)));
                if (next_bit) {
                    cur_node_rank = rank1;
                    cur_node_bv_size = cur_node_one_count;
                    ++left_siblings_count;
                    left_siblings_total_bv_size += cur_node_zero_count;
                } else {
                    cur_node_rank = rank0;
                    cur_node_bv_size = cur_node_zero_count;
                }

                // Update navigational info.
                if (depth + 1 != codelen) {
                    // Decode leaf count and total size of bitvectors corresponding to internal
                    // nodes on the next level of the tree from the (variable size) block header.
                    std::uint64_t next_level_leaf_count = *variable_block_header_ptr; ++variable_block_header_ptr;
                    std::uint64_t next_level_total_bv_size = (std::uint64_t)(*((std::uint16_t*)variable_block_header_ptr)) + 1; variable_block_header_ptr += 2;

                    // Update the total size of bitvectors corresponsding to left siblibgs of current node.
                    left_siblings_total_bv_size -= (cur_depth_total_bv_size - next_level_total_bv_size);

                    // Update bitvector offset and total size of bitvectors at current depth.
                    bv_offset += cur_depth_total_bv_size;
                    cur_depth_total_bv_size = next_level_total_bv_size;

                    // Update the number of internal nodes at current level.
                    int_nodes_count <<= 1;
                    int_nodes_count -= next_level_leaf_count;

                    // Update the number of left siblings of current node.
                    left_siblings_count -= next_level_leaf_count;
                }
            }
#else
            std::uint64_t cur_depth_bv_offset = m_superblock_header[superblock_id].m_block_headers[block_id].m_bv_offset;
            std::uint64_t cur_depth_bv_total_size = this_block_size;  // total size of bitvectors at current depth
            std::uint64_t cur_node_bv_size = this_block_size;         // size of bitvector of the current node
            std::uint64_t cur_node_rank = block_i;                    // the rank value that we refine at each level
            std::uint64_t prec_bv_total_size = 0;                     // total size of bitvectors to the left of current node at current depth
            bool cached_query = false;
            std::uint64_t cached_query_argument = 0;
            std::uint64_t cached_query_result = 0;
            for (std::uint64_t depth = 0; depth < codelen; ++depth) {
                // Issue rank query at the beginnig, end, and position `block_rank' of the current bitvector.
                // A small optimization is that sometimes the last query issued at a given depth is the same
                // as the first query on the next depth. Thus we cache the last result of rank query.
                std::uint64_t rank_beg = 0;
                if (cached_query && cached_query_argument == cur_depth_bv_offset + prec_bv_total_size)
                    rank_beg = cached_query_result;
                else rank_beg = m_superblock_header[superblock_id].m_rank_support.rank(cur_depth_bv_offset + prec_bv_total_size);
                std::uint64_t rank_mid = m_superblock_header[superblock_id].m_rank_support.rank(cur_depth_bv_offset + prec_bv_total_size + cur_node_rank);
                std::uint64_t rank_end = m_superblock_header[superblock_id].m_rank_support.rank(cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size);
                cached_query_argument = cur_depth_bv_offset + prec_bv_total_size + cur_node_bv_size;
                cached_query_result = rank_end;
                cached_query = true;

                // The number of 0/1 bits in the current internal node.
                std::uint64_t cur_node_bv_one_count = rank_end - rank_beg;
                std::uint64_t cur_node_bv_zero_count = cur_node_bv_size - cur_node_bv_one_count;

                // The number of 0/1 bits up to queried position.
                std::uint64_t rank1 = rank_mid - rank_beg;
                std::uint64_t rank0 = cur_node_rank - rank1;

                // Update rank.
                if (code & (1UL << (codelen - depth - 1))) {
                    cur_node_rank = rank1;
                    cur_node_bv_size = cur_node_bv_one_count;
                    prec_bv_total_size += cur_node_bv_zero_count;
                } else {
                    cur_node_rank = rank0;
                    cur_node_bv_size = cur_node_bv_zero_count;
                }

                // Update navigational info.
                if (depth + 1 != codelen) {
                    ++variable_block_header_ptr;
                    std::uint64_t next_level_bv_total_size = (std::uint64_t)(*((std::uint16_t*)variable_block_header_ptr)) + 1;
                    variable_block_header_ptr += 2;
                    prec_bv_total_size -= (cur_depth_bv_total_size - next_level_bv_total_size);
                    cur_depth_bv_offset += cur_depth_bv_total_size;
                    cur_depth_bv_total_size = next_level_bv_total_size;
                }
            }
#endif

            return rank_at_hyperblock_boundary + rank_at_superblock_boundary +
                   rank_at_block_boundary + cur_node_rank;
        }

        // Swap method
        void swap(wt_fbb& tree)
        {
            if (this != &tree) {
                std::swap(m_size, tree.m_size);
                std::swap(m_hyperblock_rank, tree.m_hyperblock_rank);
                std::swap(m_superblock_rank, tree.m_superblock_rank);
                std::swap(m_global_mapping, tree.m_global_mapping);
                std::swap(m_superblock_header, tree.m_superblock_header);
                std::swap(m_count, tree.m_count);
            }
        }

        // Assignment operator
        wt_fbb& operator=(const wt_fbb& tree)
        {
            if (this != &tree)
                copy(tree);
            return *this;
        }

        // Move assignment operator
        wt_fbb& operator=(wt_fbb&& tree)
        {
            swap(tree);
            return *this;
        }

        // Returns the size of the original sequence
        std::uint64_t size() const
        {
            return m_size;
        }

        // Serializes the data structure into a stream
        std::uint64_t serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            std::uint64_t written_bytes = 0;
            written_bytes += sdsl::serialize(m_size, out, child, "size");
            written_bytes += sdsl::serialize(m_count, out, child, "count");
            written_bytes += sdsl::serialize(m_hyperblock_rank, out, child, "hyperblock_rank");
            written_bytes += sdsl::serialize(m_superblock_rank, out, child, "superblock_rank");
            written_bytes += sdsl::serialize(m_global_mapping, out, child, "global_mapping");
            written_bytes += sdsl::serialize(m_superblock_header, out, child, "superblock_header");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        // Load the data structure from a stream and set the supported vector
        void load(std::istream& in)
        {
            sdsl::load(m_size, in);
            sdsl::load(m_count, in);
            sdsl::load(m_hyperblock_rank, in);
            sdsl::load(m_superblock_rank, in);
            sdsl::load(m_global_mapping, in);
            sdsl::load(m_superblock_header, in);
            for (std::uint64_t i = 0; i < m_superblock_header.size(); ++i)
                m_superblock_header[i].m_rank_support.set_vector(&m_superblock_header[i].m_bitvector);
        }
};

#ifdef ALLOW_VARIABLE_BLOCK_SIZE
template<class t_bitvector, class t_rank, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_sbs_log>::k_superblock_size_log = t_sbs_log;
template<class t_bitvector, class t_rank,std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_sbs_log>::k_superblock_size = (1UL << t_sbs_log);
template<class t_bitvector, class t_rank, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_sbs_log>::k_hyperblock_size = (1UL << 32);
#else
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_block_size_log = t_bs_log;
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_block_size = (1UL << t_bs_log);
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_superblock_size_log = t_sbs_log;
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_superblock_size = (1UL << t_sbs_log);
template<class t_bitvector, class t_rank, std::uint64_t t_bs_log, std::uint64_t t_sbs_log>
const std::uint64_t wt_fbb<t_bitvector, t_rank, t_bs_log, t_sbs_log>::k_hyperblock_size = (1UL << 32);
#endif

#endif  // __FBB_WAVELET_TREE_H_INCLUDED
