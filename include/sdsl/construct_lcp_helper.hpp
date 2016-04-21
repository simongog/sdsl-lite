#ifndef INCLUDED_SDSL_CONSTRUCT_LCP_HELPER
#define INCLUDED_SDSL_CONSTRUCT_LCP_HELPER

#include "sdsl/int_vector.hpp"
#include <queue>
#include <list>
#include <vector>

namespace sdsl
{


void insert_lcp_values(int_vector<>& partial_lcp, bit_vector& index_done, std::string lcp_file, uint64_t max_lcp_value, uint64_t lcp_value_offset);

template<class tWT>
void create_C_array(std::vector<uint64_t>& C, const tWT& wt)
{
    uint64_t quantity;                          // quantity of characters in interval
    std::vector<unsigned char> cs(wt.sigma);      // list of characters in the interval
    std::vector<uint64_t> rank_c_i(wt.sigma);    // number of occurrence of character in [0 .. i-1]
    std::vector<uint64_t> rank_c_j(wt.sigma);    // number of occurrence of character in [0 .. j-1]

    C = std::vector<uint64_t>(257, 0);
    interval_symbols(wt, 0, wt.size(), quantity, cs, rank_c_i, rank_c_j);
    for (uint64_t i=0; i<quantity; ++i) {
        unsigned char c = cs[i];
        C[c+1] = rank_c_j[i];
    }
    for (uint64_t i=1; i<C.size()-1; ++i) {
        C[i+1] += C[i];
    }
}


class buffered_char_queue
{
        typedef bit_vector::size_type size_type;
        typedef std::queue<uint8_t> tQ;
    private:
        static const uint32_t m_buffer_size =  10000;//409600;
        uint8_t m_write_buf[m_buffer_size];
        uint8_t m_read_buf[m_buffer_size];
        size_type 	m_widx; // write index
        size_type 	m_ridx; // read index
        bool		m_sync; // are read and write buffer the same?
        size_type 	m_disk_buffered_blocks; // number of blocks written to disk and not read again yet
        char 		m_c;
        size_type	m_rb; // read blocks
        size_type	m_wb; // written blocks

        std::string m_file_name;

        std::fstream	m_stream;

    public:

        buffered_char_queue();
        void init(const std::string& dir, char c);
        ~buffered_char_queue();
        void push_back(uint8_t x);
        uint8_t pop_front();
};

typedef std::list<int_vector<>::size_type> tLI;
typedef std::vector<int_vector<>::size_type> tVI;

template<class size_type_class>
void push_front_m_index(size_type_class i, uint8_t c, tLI(&m_list)[256], uint8_t (&m_chars)[256], size_type_class& m_char_count)
{
    if (m_list[c].empty()) {
        m_chars[m_char_count++] = c;
    }
    m_list[c].push_front(i);
}

template<class size_type_class>
void push_back_m_index(size_type_class i, uint8_t c, tLI(&m_list)[256], uint8_t (&m_chars)[256], size_type_class& m_char_count)
{
    if (m_list[c].empty()) {
        m_chars[m_char_count++] = c;
    }
    m_list[c].push_back(i);
}

void lcp_info(tMSS& file_map);

}

#endif
