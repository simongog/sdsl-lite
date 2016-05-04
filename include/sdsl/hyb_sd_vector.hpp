/* sdsl - succinct data structures library
    Copyright (C) 2016 Simon Gog, Matthias Petri

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
/*!\file sd_vector.hpp
   \brief sd_vector.hpp contains the sdsl::sd_vector class, and
          classes which support rank and select for sd_vector.
   \author Simon Gog, Matthias Petri
*/
#ifndef INCLUDED_SDSL_HYB_SD_VECTOR
#define INCLUDED_SDSL_HYB_SD_VECTOR

#include "int_vector.hpp"
#include "sd_vector.hpp"
#include "util.hpp"
#include "iterators.hpp"

//! Namespace for the succinct data structure library
namespace sdsl
{


constexpr bool debug=false;
size_t g_hi_size=0;

template <class t_itr>
std::string print_vec(t_itr beg, t_itr end)
{
    std::string str = "[";
    auto itr = beg;
    while (itr != (end - 1)) {
        str += std::to_string(*itr) + " ";
        ++itr;
    }
    str += std::to_string(*itr) + "]";
    return str;
}

inline uint64_t next0(const uint64_t* word, uint64_t idx)
{
    word += (idx >> 6);
    auto masked_inverse_word = ~(*word | sdsl::bits::lo_set[(idx & 0x3F) + 1]);
    if (masked_inverse_word) {
        return (idx & ~((size_t)0x3F)) + sdsl::bits::lo(masked_inverse_word);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    ++word;
    while (*word == 0xFFFFFFFFFFFFFFFFULL) {
        idx += 64;
        ++word;
    }
    return idx + sdsl::bits::lo(~(*word));
}

/*!
 * \param word Beginning of bit_vector (represented as sequence of uint64_t words)
 * \param idx  Initial scanning position (in bits)
 * \param i    i
 * \return The number of set bits up to position i (exlusive)
 */
template <uint16_t t_block_size>
inline uint64_t cnt(const uint64_t* word, uint64_t idx, uint64_t i)
{
//    std::cout<<"cnt("<<idx<<", "<<i<<")"<<std::endl;
    word += (idx >> 6);
    auto offset = idx & 0x3F;
//    std::cout<<"offset="<<offset<<std::endl;
    uint64_t w = (*word) >> offset;
    uint64_t pre_considered = 0;
    uint64_t considered = 64 - offset;
//    std::cout<<"considered="<<considered<<std::endl;
    uint64_t res = 0;
    while (considered < i) {
        res += sdsl::bits::cnt(w);
//        std::cout<<"res="<<res<<" after "<<considered<<" bits"<<std::endl;
        if (res >= t_block_size) {
            return t_block_size;
        }
        pre_considered = considered;
        considered += 64;
        w = *(++word);
    }

//std::cout<<"considered="<<considered<<std::endl;
    // considered \in [i+0..i+63]
    if (i == considered) {
        return res + sdsl::bits::cnt(w);
    }
//std::cout<<"pre_considered="<<pre_considered<<std::endl;
    i = i-pre_considered;//i + 64 - considered;
    w &= sdsl::bits::lo_set[i];
    res += sdsl::bits::cnt(w);
    return res > t_block_size ? t_block_size : res; // TODO: is the space reduction large enough to justify this  ???
}
/*!
 * \param word Beginning of bit_vector (represented as sequence of uint64_t words)
 * \param idx  Initial scanning position (in bits)
 * \param i    i
 * \return The absolut position (in bits) of the i-th set bit (\f$ i>0 \f$ from idx
 */
inline uint64_t sel(const uint64_t* word, uint64_t idx, uint64_t i)
{
    --i;
    word += (idx >> 6);
    auto masked_word = *word & ~sdsl::bits::lo_set[(idx & 0x3F) + 1];
    auto one_cnt = sdsl::bits::cnt(masked_word);
    if (one_cnt >= i) {
        return (idx & ~((size_t)0x3F)) + sdsl::bits::sel(masked_word, i);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    i -= one_cnt;
    ++word;
    one_cnt = sdsl::bits::cnt(*word);
    while (i > one_cnt) {
        ++word;
        idx += 64;
        i -= one_cnt;
        one_cnt = sdsl::bits::cnt(*word);
    }
    return idx + sdsl::bits::sel(*word, i);
}

/*!
 * \param word Beginning of bit_vector (represented as sequence of uint64_t words)
 * \param idx  Initial scanning position (in bits)
 * \param i    i
 * \return The absolut position (in bits) of the i-th unset bit (\f$ i>0 \f$ from idx
 */
inline uint64_t sel0(const uint64_t* word, uint64_t idx, uint64_t i, uint64_t max_considered)
{
    //    std::cout<<"cnt("<<idx<<", "<<i<<")"<<std::endl;
    word += (idx >> 6);
    auto offset = idx & 0x3F;
    uint64_t w = (~(*word)) >> offset;
    uint64_t considered = 64 - offset;
    uint64_t res = 0;
    uint64_t cnt = 0;
    uint64_t word_cnt = sdsl::bits::cnt(w);

    while (cnt + word_cnt < i) {
        cnt += word_cnt;
        if (considered > max_considered) {
//             std::cout << "LOOP considered = " << considered << " max_considered = " << max_considered << std::endl;
            return std::numeric_limits<uint64_t>::max();
        }
        res = considered;
        considered += 64;
        w = (~(*(++word)));
        word_cnt = sdsl::bits::cnt(w);
    }
    // cnt < i and cnt+word_cnt >= i
    // add select (i-cnt) to res
    res += sdsl::bits::sel(w, i - cnt);
    if (res  > max_considered) {
//        std::cout << "res = " << res << " i = " << i << std::endl;
//        std::cout << "IF considered = " << considered << " max_considered = " << max_considered << std::endl;
        return std::numeric_limits<uint64_t>::max();
    }
    return res;
}

template <uint64_t t_block_size>
class hyb_sd_block_support_bv
{
    public:
        typedef bit_vector::size_type size_type;
        typedef size_type value_type;

    public:
        static size_type estimate_size(size_type u)
        {
            return u;
        }
        static size_type
        serialize(sdsl::bit_vector& bv, size_type offset, sdsl::int_vector<64>& data, size_type)
        {
//        std::cout<<"write BV of size "<<data[t_block_size-1]+1<<std::endl;
            for (size_t i = 0; i < data[t_block_size - 1]+1; ++i)
                bv[offset + i] = 0;
            for (size_type i = 0; i < t_block_size; ++i) {
                bv[offset + data[i]] = 1;
            }
//        std::cout<<"serialize bv: "<<data[t_block_size-1]+1<<std::endl;
//        std::cout<<"bv="<<bv<<std::endl;
//        std::cout<<"data="<<data<<std::endl;
            return data[t_block_size - 1] + 1;
        }

        static size_type select_1(const sdsl::bit_vector& bv, size_type offset, size_type i, size_type)
        {
            return sel(bv.data(), offset, i + 1) - offset;
        }

        static size_type rank_1(const sdsl::bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type)
        {
            auto offset = block_start[block_id];
            auto next_offset = block_start[block_id+1];
            if (i > next_offset-offset)
                return t_block_size;
            /*
                    std::cout<<"offset="<<offset<<" i="<<i<<" block_id="<<block_id<<std::endl;
                    size_t myrank=0;
                    for(size_t j=offset; j<offset+i; ++j){
                        myrank += bv[j];
                    }
                    std::cout<<"myrank="<<myrank<<std::endl;
                    size_t rank = cnt<t_block_size>(bv.data(), offset, i);
                    std::cout<<"rank="<<rank<<std::endl;
            */
            return cnt<t_block_size>(bv.data(), offset, i);
        }
};

template <uint64_t t_block_size>
class hyb_sd_block_support_full
{
    public:
        typedef bit_vector::size_type size_type;
        typedef size_type value_type;

    public:
        static size_type estimate_size(size_type u)
        {
            return (t_block_size == u) ? 0 : std::numeric_limits<size_type>::max();
        }

        static size_type
        serialize(sdsl::bit_vector&, size_type, sdsl::int_vector<64>&, size_type u)
        {
            if (t_block_size != u) {
                std::cerr << "this should not happen!" << std::endl;
            }
            return 0;
        }

        static size_type select_1(const sdsl::bit_vector&, size_type, size_type i, size_type)
        {
            return i;
        }

        static size_type rank_1(const sdsl::bit_vector&, const int_vector<>&, size_type, size_type i, size_type)
        {
            return i;
        }
};

template <uint64_t t_block_size>
class hyb_sd_block_support_ef
{
    public:
        typedef bit_vector::size_type size_type;
        typedef size_type value_type;

    private:
        // TODO factor out calculation of logu and logm

    public:
        static size_type estimate_size(size_type u)
        {
            uint8_t logu = sdsl::bits::hi(u) + 1;
            uint8_t logm = sdsl::bits::hi(t_block_size) + 1; // TODO constexpr for hi?
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;
            size_type size_in_bits = width_low * t_block_size + (1ULL << logm) + t_block_size + 1;
            return size_in_bits;
        }
        static size_type
        serialize(sdsl::bit_vector& bv, size_type offset, sdsl::int_vector<64>& data, size_type u)
        {
            size_type written_bits = 0;
            uint8_t logu = sdsl::bits::hi(u) + 1;
            uint8_t logm = sdsl::bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;
//std::cout<<"offst="<<offset<<" u="<<u<<" data="<<data<<" width_low="<<width_low<<std::endl;
            /* write low */
            auto data_ptr = bv.data() + (offset / 64);
            uint8_t in_word_offset = offset % 64;
            for (size_type i = 0; i < t_block_size; i++) {
                uint64_t x = data[i];
                sdsl::bits::write_int_and_move(data_ptr, x, in_word_offset, width_low);
            }
            written_bits += width_low * t_block_size;

            /* write high */
            size_type last_high = 0;
            for (size_type i = 0; i < t_block_size; i++) {
                uint64_t x = data[i];
                size_type cur_high = x >> width_low;
                size_type write_val = cur_high - last_high;
                while (write_val >= 64) {
                    sdsl::bits::write_int_and_move(data_ptr, 0ULL, in_word_offset, 64);
                    write_val -= 64;
                    written_bits += 64;
                }
                sdsl::bits::write_int_and_move(data_ptr, 1ULL << write_val, in_word_offset, write_val + 1);
                last_high = cur_high;
                written_bits += write_val + 1;
            }
            bv[offset+written_bits] = 0;
            ++written_bits;
            return written_bits;
        }

        static size_type select_1(const sdsl::bit_vector& bv, size_type offset, size_type i, size_type u)
        {
            uint8_t logu = sdsl::bits::hi(u) + 1;
            uint8_t logm = sdsl::bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;
            size_type hi_part_offset = offset + t_block_size * width_low;

            if (debug) {
                std::cout<<"i="<<i<<" u="<<u<<" logu="<<(int)logu<<" width_low="<<width_low<<" logm="<<(int)logm<<std::endl;

                size_t num_one = 0;
                size_t num_zero = 0;
                for (size_t i=0; i<g_hi_size; i++) {
                    std::cout << "hi["<<i<<"] = " << bv[hi_part_offset+i];
                    if (bv[hi_part_offset+i]) {
                        std::cout << "  one # " << ++num_one;
                    } else {
                        std::cout << " zero # " << ++num_zero;
                    }
                    std::cout << std::endl;
                }
                std::cout << "hi_size = " << g_hi_size<<std::endl;
            }
            size_type low_part_offset = offset + i * width_low;

            auto low_part_data_ptr = bv.data() + (low_part_offset / 64);
            uint8_t low_part_in_word_offset = low_part_offset % 64;
            auto low_part = sdsl::bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low);

            auto bucket = sel(bv.data(), hi_part_offset, i + 1) - hi_part_offset - i;
            if (debug) {
                std::cout<<"bucket="<<bucket<<" sel(bv.data(),"<<hi_part_offset<<","<<i+1<<")="<<sel(bv.data(),hi_part_offset,i+1) <<std::endl;
            }
            return (bucket << width_low) | low_part;
        }

        static size_type rank_1(const sdsl::bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type u)
        {
            auto offset = block_start[block_id];
            auto next_offset = block_start[block_id + 1];

            uint8_t logu = sdsl::bits::hi(u) + 1;
            uint8_t logm = sdsl::bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;

            size_type hi_part_offset = offset + t_block_size * width_low;
            size_type hi_size = next_offset - hi_part_offset;

//std::cout<<"i="<<i<<" u="<<u<<" logu="<<(int)logu<<" width_low="<<width_low<<" logm="<<(int)logm<<std::endl;
//
//size_t num_one = 0;
// size_t num_zero = 0;
// for(size_t i=0;i<hi_size;i++) {
//     std::cout << "hi["<<i<<"] = " << bv[hi_part_offset+i];
//     if(bv[hi_part_offset+i]) {
//         std::cout << "  one # " << ++num_one;
//     } else {
//         std::cout << " zero # " << ++num_zero;
//     }
//     std::cout << std::endl;
// }
// std::cout << "hi_size = " << hi_size<<std::endl;

            size_type high_val = (i >> width_low);
            size_type local_sel = sel0(bv.data(), hi_part_offset, high_val + 1, hi_size);
            if (local_sel == std::numeric_limits<uint64_t>::max()) {
//            std::cout<<high_val<<std::endl;
//            std::cout<<"local_sel="<<local_sel<<std::endl;
                return t_block_size;
            }

            size_type sel_high = local_sel;
            size_type rank_low = sel_high - high_val;
            if (0 == rank_low) {
//            std::cout<<"bye"<<std::endl;
                return 0;
            }

            size_type low_part_offset = offset + rank_low * width_low;
            size_type val_low = i & bits::lo_set[width_low];
//std::cout << "RANK1-EF(i=" << i << ",u=" << u << ",wl=" << width_low << ",sel_high=" << sel_high << ",high_val=" << high_val << ",rank_low=" << rank_low << ",val_low=" << val_low << ")" << std::endl;

            auto low_part_data_ptr = bv.data() + (low_part_offset / 64);
            uint8_t low_part_in_word_offset = low_part_offset % 64;

            do {
                if (!sel_high)
                    return 0;
                --rank_low;
                --sel_high;
                low_part_offset -= width_low;
                low_part_data_ptr = bv.data() + (low_part_offset / 64);
                low_part_in_word_offset = low_part_offset % 64;
//std::cout << "cur_low ("<<rank_low<<") = " << sdsl::bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low) << std::endl;
            } while (bv[hi_part_offset + sel_high] and sdsl::bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low) >= val_low);

//std::cout << "AFER SCAN RANK1-EF(i=" << i << ",u=" << u << ",wl=" << width_low << ",sel_high=" << sel_high << ",high_val=" << high_val << ",rank_low=" << rank_low << ",val_low=" << val_low << ")" << std::endl;

            return rank_low + 1;
        }
};

template <uint8_t t_b,
          class t_hyb_sd_bv>
class select_support_hyb_sd;

template <uint8_t t_b,
          class t_hyb_sd_bv>
class rank_support_hyb_sd;

template <uint16_t t_block_size = 128>
class hyb_sd_vector
{
    public:
        typedef bit_vector::size_type size_type;
        typedef bool value_type;
        typedef bit_vector::difference_type difference_type;
        typedef random_access_const_iterator<hyb_sd_vector> iterator;
        typedef iterator const_iterator;
        typedef bv_tag index_category;

        // typedef rank_support_hyb_sd<0, hyb_sd_vector> rank_0_type;
        typedef rank_support_hyb_sd<1, hyb_sd_vector> rank_1_type;
        // typedef select_support_hyb_sd<0, hyb_sd_vector> select_0_type;
        typedef select_support_hyb_sd<1, hyb_sd_vector> select_1_type;

    private:
        sd_vector<> m_top;
        sd_vector<>::select_1_type m_top_sel;
        sd_vector<>::rank_1_type m_top_rank;
        sdsl::bit_vector m_bottom;
        sdsl::int_vector<> m_block_start;
        size_type m_size = 0;
        size_type m_num_ones = 0;


    public:
        static constexpr uint16_t block_size = t_block_size;

    private:
        enum class hyb_sd_blocktype
        {
            EF,
            BV,
            FULL
        };

        std::pair<hyb_sd_blocktype, size_type>
        determine_block_type(size_t u) const
        {
            if (u == t_block_size)
                return { hyb_sd_blocktype::FULL, 0 };
            size_type ef_bits = hyb_sd_block_support_ef<t_block_size>::estimate_size(u);
            size_type bv_bits = hyb_sd_block_support_bv<t_block_size>::estimate_size(u);
            if (bv_bits < ef_bits)
                return { hyb_sd_blocktype::BV, bv_bits };
            return { hyb_sd_blocktype::EF, ef_bits };
        }

        size_type compress_block(size_type offset, sdsl::int_vector<64>& data, size_t u)
        {
            size_type written_bits = 0;
            auto bt = determine_block_type(u);
            auto type = bt.first;
            auto size_in_bits = bt.second;
            if (m_bottom.size() < offset + size_in_bits) {
                m_bottom.resize(m_bottom.size() * 2 + size_in_bits);
            }

            for (size_t i=1; i<data.size(); i++) {
                if (data[i] == data[i-1]) {
                    std::cout << "error: equal data encoded i="<<i<< " d[i]="<<data[i]<< " d[i-1]="<<data[i-1] << std::endl;
                }
            }

            switch (type) {
                case hyb_sd_blocktype::BV:
                    written_bits = hyb_sd_block_support_bv<t_block_size>::serialize(m_bottom, offset, data, u);
                    break;
                case hyb_sd_blocktype::EF:
                    written_bits = hyb_sd_block_support_ef<t_block_size>::serialize(m_bottom, offset, data, u);
                    break;
                case hyb_sd_blocktype::FULL:
                    /* nothing to store */
                    break;
            }
            return written_bits;
        }

    public:

        hyb_sd_vector() {}

///*
        explicit hyb_sd_vector(const bit_vector& bv) //: hyb_sd_vector(bv.ones_begin(),bv.ones_end(),bv.size())
        {
            if (bv.size() == 0) {
                return;
            }
            m_size = bv.size();
            m_num_ones = select_support_trait<1, 1>::arg_cnt(bv);
            if (m_num_ones == 0) {
                return;
            }
            size_type num_full_blocks = m_num_ones / t_block_size;
            size_type num_blocks = num_full_blocks;
            size_type num_leftover = m_num_ones % t_block_size;
            size_type num_dummy = (t_block_size - num_leftover) % t_block_size;
            bool has_leftover_block = num_leftover != 0;
            if (has_leftover_block) {
                num_blocks++;
            }

            // (1) fill the top level
            std::vector<uint64_t> top_lvl;
            size_t one_found = t_block_size - 1; // we always want to add the first one!
            size_type last_one = 0;
            for (size_type i = 0; i < bv.size(); i++) {
                if (bv[i] == 1) {
                    last_one = i;
                    if (++one_found == t_block_size) {
                        top_lvl.push_back(i);
                        one_found = 0;
                    }
                }
            }
            // terminate the top level so top[i+1] - top[i] always works
//        top_lvl.push_back(last_one + num_dummy + 1);
            top_lvl.push_back(std::max(bv.size(), last_one+1) + num_dummy);

            // (2) bottom level
            m_block_start.resize(num_blocks + 1);
            size_type value_offset = 0;
            size_type written_bits = 0;
            sdsl::int_vector<64> tmp_data(t_block_size);
            size_t j = 0;
            for (size_type i = 0; i < num_blocks; i++) {
                m_block_start[i] = written_bits;
                // (2a) compute block data
                value_offset = top_lvl[i];
                one_found = 0;
                while (one_found != t_block_size) {
                    if (j <= last_one) {
                        if (bv[j] == 1) {
                            tmp_data[one_found++] = j - value_offset;
                        }
                    }
//                else if ( j == last_one+1 ) {
//                    tmp_data[one_found] = m_size - value_offset;
//                    one_found++;
//                }
                    else {
                        tmp_data[one_found] = std::max(m_size-value_offset, tmp_data[one_found - 1] + 1);
//                    std::cout<<"pad "<<tmp_data[one_found]+value_offset<<std::endl;
                        one_found++;
                    }
                    ++j;
                }
                // (2b) compress block
                size_type block_universe = top_lvl[i + 1] - top_lvl[i];
                auto wb = compress_block(m_block_start[i], tmp_data, block_universe);
                written_bits += wb;
            }
            m_block_start[num_blocks] = written_bits;
            m_bottom.resize(written_bits);

            // (3) encode the top level
            m_top = decltype(m_top)(top_lvl.begin(), top_lvl.end());
            m_top_sel = decltype(m_top_sel)(&m_top);
            m_top_rank = decltype(m_top_rank)(&m_top);

            // (4) bit compress pointers
            sdsl::util::bit_compress(m_block_start);
        }
//*/

        template <class t_itr>
        hyb_sd_vector(const t_itr begin, const t_itr end, size_type bv_size = 0)
        {
            if (begin == end and bv_size==0) {
                return;
            }
            if (!is_sorted(begin, end)) {
                throw std::runtime_error("hyb_sd_vector: source list is not sorted.");
            }
            m_size = bv_size;
            if (bv_size == 0)
                m_size = *(end - 1) + 1;
            m_num_ones = std::distance(begin, end);
            if (m_num_ones==0) {
                return;
            }
//std::cout<<"m_num_ones="<<m_num_ones<<std::endl;
            size_type num_full_blocks = m_num_ones / t_block_size;
            size_type num_blocks = num_full_blocks;
            size_type num_leftover = m_num_ones % t_block_size;
            size_type num_dummy = (t_block_size - num_leftover) % t_block_size;
            bool has_leftover_block = num_leftover != 0;
            if (has_leftover_block) {
                num_blocks++;
            }

            // (1) fill the top level
            std::vector<uint64_t> top_lvl;
            auto itr = begin;
            while (itr < end) {
                top_lvl.push_back(*itr);
                itr += t_block_size;
            }
            // terminate the top level so top[i+1] - top[i] always works
            top_lvl.push_back(std::max(bv_size, *(end - 1)+1) + num_dummy);

            // (2) bottom level
            m_block_start.resize(num_blocks + 1);
//std::cout<<"num_blocks="<<num_blocks<<std::endl;
            itr = begin;
            size_type value_offset = 0;
            size_type written_bits = 0;
            sdsl::int_vector<64> tmp_data(t_block_size);
            for (size_type i = 0; i < num_blocks; i++) {
                m_block_start[i] = written_bits;
                // (2a) compute block data
                value_offset = top_lvl[i];
                for (size_type j = 0; j < t_block_size; j++) {
                    if (itr == end) {
                        tmp_data[j] = std::max(m_size-value_offset, tmp_data[j - 1] + 1);
                    } else {
                        tmp_data[j] = *itr - value_offset;
                        ++itr;
                    }
                }

                // (2b) compress block
                size_type block_universe = top_lvl[i + 1] - top_lvl[i];
                // std::cout << "compress block " << i << std::endl;
                auto wb = compress_block(m_block_start[i], tmp_data, block_universe);
                written_bits += wb;
            }
            m_block_start[num_blocks] = written_bits;
            m_bottom.resize(written_bits);

            // (3) encode the top level
            m_top = decltype(m_top)(top_lvl.begin(), top_lvl.end());
            m_top_sel = decltype(m_top_sel)(&m_top);
            m_top_rank = decltype(m_top_rank)(&m_top);

            // (4) bit compress pointers
            sdsl::util::bit_compress(m_block_start);
//std::cout<<"bye hyb_sd"<<std::endl;
        }


        value_type operator[](size_type i)const
        {
            return rank_1(i+1) - rank_1(i);
        }

        //! Accessing the i-th element of the original bit_vector
        size_type select_1(size_type i) const
        {
//debug = 23063 == i;
            i = i - 1;
            auto block_id = i / t_block_size;
            auto in_block_offset = i % t_block_size;
            auto top_value = m_top_sel(block_id + 1);
            size_type res = top_value;

            if (in_block_offset == 0)
                return res;

            auto u = m_top_sel(block_id + 2) - top_value;
            auto bt = determine_block_type(u);
            auto block_type = bt.first;
            size_type block_offset = m_block_start[block_id];

            if (debug) {
                std::cout<<"block_id = " << block_id << std::endl;
                std::cout<<"block_offset = " << block_offset << std::endl;
                std::cout<<"in_block_offset = " << in_block_offset << std::endl;
                std::cout<<"res = " << res << std::endl;
                std::cout<<"u = " << u << std::endl;
            }

            switch (block_type) {
                case hyb_sd_blocktype::BV:
                    if (debug) std::cout << "BV" << std::endl;
                    res += hyb_sd_block_support_bv<t_block_size>::select_1(m_bottom, block_offset, in_block_offset, u);
                    break;
                case hyb_sd_blocktype::EF:
                    if (debug) {
                        std::cout << "EF" << std::endl;
                        g_hi_size = m_block_start[block_id+1]-m_block_start[block_id];
                    }
                    res += hyb_sd_block_support_ef<t_block_size>::select_1(m_bottom, block_offset, in_block_offset, u);
                    break;
                case hyb_sd_blocktype::FULL:
                    if (debug) std::cout << "FULL" << std::endl;
                    res += hyb_sd_block_support_full<t_block_size>::select_1(m_bottom, block_offset, in_block_offset, u);
                    break;
            }
            return res;
        }

        size_type rank_1(size_type i) const
        {
            if (i > m_size or m_num_ones == 0) {
                return m_num_ones;
            }
            auto block_id = m_top_rank(i);
            if (block_id == 0)
                return 0;
            block_id -= 1;
            size_type res = block_id * t_block_size;
            auto top_value = m_top_sel(block_id + 1);
            size_type in_block_i = i;
            in_block_i -= top_value;
            if (in_block_i == 0)
                return res;

            auto u = m_top_sel(block_id + 2) - top_value;
            auto bt = determine_block_type(u);
            auto block_type = bt.first;
            size_type block_offset = m_block_start[block_id];
            switch (block_type) {
                case hyb_sd_blocktype::BV:
//             std::cout << "BV" << std::endl;
                    res += hyb_sd_block_support_bv<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
                case hyb_sd_blocktype::EF:
//             std::cout << "EF" << std::endl;
                    res += hyb_sd_block_support_ef<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
                case hyb_sd_blocktype::FULL:
//             std::cout << "FULL" << std::endl;
                    res += hyb_sd_block_support_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
            }

            return res;
        }

        //! Get the integer value of the binary string of length len starting at position idx.
        uint64_t get_int(size_type idx, const uint8_t len = 64) const
        {
            uint64_t x = 0ULL;
            for (size_t i=0; i<len and idx+i < size(); ++i) {
                x |= (static_cast<uint64_t>((*this)[idx+i])) << i;
            }
            return x;
        }

        //! Returns the size of the original bit vector.
        size_type size() const
        {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_num_ones, out, child, "num_ones");
            written_bytes += m_top.serialize(out, child, "top");
            written_bytes += m_top_sel.serialize(out, child, "top_sel");
            written_bytes += m_top_rank.serialize(out, child, "top_rank");
            written_bytes += m_bottom.serialize(out, child, "bottom");
            written_bytes += m_block_start.serialize(out, child, "block_start");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            read_member(m_num_ones, in);
            m_top.load(in);
            m_top_sel.load(in);
            m_top_sel.set_vector(&m_top);
            m_top_rank.load(in);
            m_top_rank.set_vector(&m_top);
            m_bottom.load(in);
            m_block_start.load(in);
        }

        void swap(hyb_sd_vector& v)
        {
            std::swap(m_size, v.m_size);
            std::swap(m_num_ones, v.m_num_ones);
            m_top.swap(v.m_top);
            util::swap_support(m_top_sel, v.m_top_sel, &m_top, &(v.m_top));
            util::swap_support(m_top_rank, v.m_top_rank, &m_top, &(v.m_top));
            m_bottom.swap(v.m_bottom);
            m_block_start.swap(v.m_block_start);
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }
};

//! Select data structure for hyb_sd_vector
template <uint8_t t_b = 1,
          class hyb_bv_type = hyb_sd_vector<>>
class select_support_hyb_sd
{
    public:
        typedef typename hyb_bv_type::size_type size_type;
        typedef hyb_bv_type bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
        static constexpr uint16_t block_size = hyb_bv_type::block_size;

    private:
        const hyb_bv_type* m_v;

    public:
        explicit select_support_hyb_sd(const hyb_bv_type* v = nullptr)
        {
            set_vector(v);
        }

        size_type select(size_type i) const
        {
            return m_v->select_1(i);
        }

        size_type operator()(size_type i) const
        {
            return select(i);
        }

        size_type size() const
        {
            return m_v->size();
        }

        void set_vector(const hyb_bv_type* v = nullptr)
        {
            m_v = v;
        }

        select_support_hyb_sd& operator=(const select_support_hyb_sd& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(select_support_hyb_sd&) {}

        void load(std::istream&, const hyb_bv_type* v = nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
        {
            return serialize_empty_object(out, v, name, this);
        }
};

//! Rank data structure for hyb_sd_vector
template <uint8_t t_b = 1,
          class hyb_bv_type = hyb_sd_vector<>>
class rank_support_hyb_sd
{
    public:
        typedef typename hyb_bv_type::size_type size_type;
        typedef hyb_bv_type bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
        static constexpr uint16_t block_size = hyb_bv_type::block_size;

    private:
        const hyb_bv_type* m_v;

    public:
        explicit rank_support_hyb_sd(const hyb_bv_type* v = nullptr)
        {
            set_vector(v);
        }

        size_type rank(size_type i) const
        {
            return m_v->rank_1(i);
        }

        size_type operator()(size_type i) const
        {
            return rank(i);
        }

        size_type size() const
        {
            return m_v->size();
        }

        void set_vector(const hyb_bv_type* v = nullptr)
        {
            m_v = v;
        }

        rank_support_hyb_sd& operator=(const rank_support_hyb_sd& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(rank_support_hyb_sd&) {}

        void load(std::istream&, const hyb_bv_type* v = nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
        {
            return serialize_empty_object(out, v, name, this);
        }
};

} // end namespace
#endif
