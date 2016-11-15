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
#include "coder.hpp"
#include "util.hpp"
#include "iterators.hpp"

//! Namespace for the succinct data structure library
namespace sdsl
{

//std::vector<uint64_t> g_range_stats;
extern size_t g_saved_bits;

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
    auto masked_inverse_word = ~(*word | bits::lo_set[(idx & 0x3F) + 1]);
    if (masked_inverse_word) {
        return (idx & ~((size_t)0x3F)) + bits::lo(masked_inverse_word);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    ++word;
    while (*word == 0xFFFFFFFFFFFFFFFFULL) {
        idx += 64;
        ++word;
    }
    return idx + bits::lo(~(*word));
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
        res += bits::cnt(w);
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
        return res + bits::cnt(w);
    }
//std::cout<<"pre_considered="<<pre_considered<<std::endl;
    i = i-pre_considered;//i + 64 - considered;
    w &= bits::lo_set[i];
    res += bits::cnt(w);
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
    auto masked_word = *word & ~bits::lo_set[(idx & 0x3F) + 1];
    auto one_cnt = bits::cnt(masked_word);
    if (one_cnt >= i) {
        return (idx & ~((size_t)0x3F)) + bits::sel(masked_word, i);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    i -= one_cnt;
    ++word;
    one_cnt = bits::cnt(*word);
    while (i > one_cnt) {
        ++word;
        idx += 64;
        i -= one_cnt;
        one_cnt = bits::cnt(*word);
    }
    return idx + bits::sel(*word, i);
}

/*!
 * \param word Beginning of bit_vector (represented as sequence of uint64_t words)
 * \param idx  Initial scanning position (in bits)
 * \param i    i
 * \return The absolut position (in bits) of the i-th unset bit (\f$ i>0 \f$ from idx
 */
inline uint64_t sel0(const uint64_t* word, uint64_t idx, uint64_t i)
{
    //    std::cout<<"cnt("<<idx<<", "<<i<<")"<<std::endl;
    word += (idx >> 6);
    auto offset = idx & 0x3F;
    uint64_t w = (~(*word)) >> offset;
    uint64_t considered = 64 - offset;
    uint64_t res = 0;
    uint64_t cnt = 0;
    uint64_t word_cnt = bits::cnt(w);

    while (cnt + word_cnt < i) {
        cnt += word_cnt;
        res = considered;
        considered += 64;
        w = (~(*(++word)));
        word_cnt = bits::cnt(w);
    }
    // cnt < i and cnt+word_cnt >= i
    // add select (i-cnt) to res
    res += bits::sel(w, i - cnt);
    return res;
}

template <uint64_t t_block_size>
class hyb_sd_block_bv
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
        serialize(bit_vector& bv, size_type offset, int_vector<64>& data, size_type)
        {
            for (size_t i = 0; i < data[t_block_size - 1]+1; ++i)
                bv[offset + i] = 0;
            for (size_type i = 0; i < t_block_size; ++i) {
                bv[offset + data[i]] = 1;
            }
            return data[t_block_size - 1] + 1;
        }

        static size_type select_1(const bit_vector& bv, size_type offset, size_type i, size_type)
        {
            return sel(bv.data(), offset, i + 1) - offset;
        }

        static size_type
        rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type)
        {
            auto offset = block_start[block_id];
            auto next_offset = block_start[block_id+1];
            if (i > next_offset-offset)
                return t_block_size;
            return cnt<t_block_size>(bv.data(), offset, i);
        }

        static std::array<size_type,2>
        rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, const std::array<size_type,2>& ij, size_type)
        {
            auto offset = block_start[block_id];
            auto next_offset = block_start[block_id+1];
            if (ij[0] > next_offset-offset) {
                return {t_block_size,t_block_size};
            }
            auto resi =  cnt<t_block_size>(bv.data(), offset, ij[0]);
            if (ij[1] > next_offset-offset) {
                return {resi, t_block_size};
            }
            return {resi, resi+cnt<t_block_size>(bv.data(), offset+ij[0], ij[1]-ij[0])};
        }

};


template <uint64_t t_block_size, typename t_coder=coder::elias_delta>
class hyb_sd_block_rl
{
    public:
        typedef bit_vector::size_type size_type;
        typedef size_type value_type;

    public:
        static size_type encode(int_vector<64>& data, bit_vector* bv=nullptr, size_type offset=0)
        {
            uint64_t* data_ptr = nullptr;
            uint8_t in_word_offset = offset % 64;
            if (bv != nullptr) {
                data_ptr = bv->data() + (offset / 64);
            }
            auto do_encode = [&](uint64_t x) {
                if (data_ptr != nullptr) {
                    t_coder::encode(x, data_ptr, in_word_offset);
                }
                return t_coder::encoding_length(x);
            };

            size_type rl_bits = 0;
            size_t begin = 0, end = 1;
            while (end < data.size()) {
                uint64_t delta = data[end]-data[begin];
                if (delta > end-begin) {
                    if (end-begin == 1) {
                        rl_bits += do_encode(delta);
                    } else { // end-begin > 1
                        rl_bits += do_encode(1);
                        rl_bits += do_encode(end-1-begin);
                        rl_bits += do_encode(data[end]-data[end-1]);
                    }
                    begin = end;
                    ++end;
                } else {
                    ++end;
                }
            }
            if (end-begin > 1) {
//                rl_bits += do_encode(1);
//                rl_bits += do_encode(end-1-begin);
                if (bv != nullptr) {
                    g_saved_bits += t_coder::encoding_length(1);
                    g_saved_bits += t_coder::encoding_length(end-1-begin);
                }
            }
            /*
                        if ( bv!=nullptr ) {
            //                std::cout<<"Checking block "<<std::endl;
            //                std::cout<<"rl_bits="<<rl_bits<<std::endl;
            //                std::cout<<"data="<<data<<std::endl;
                            in_word_offset = offset % 64;
                            data_ptr = bv->data() + (offset / 64);
                            auto end_offset = offset + rl_bits;
                            auto temp = decode(data_ptr, in_word_offset, bv->data() + (end_offset/64), end_offset%64);
                            if ( data.size() != temp.size() ){
                                std::cout<<"Error in RL block; size is different!!!"<<std::endl;
                                throw std::logic_error("error in RL block");
                            }
                            for(size_t i=0; i<data.size(); ++i){
                                if(data[i]!=temp[i]) {
                                    std::cout<<"Error in RL block; decoded value is different!!!"<<std::endl;
                                    std::cout<<"data["<<i<<"]="<< data[i] <<" != "<<temp[i]<<" = temp[i]"<<std::endl;
                                    throw std::logic_error("error in RL block");
                                }
                            }
                        }
            */
            return rl_bits;
        }

        static int_vector<64> decode(const uint64_t* data_ptr, uint8_t offset, const uint64_t* data_ptr_end, uint8_t offset_end)
        {
            int_vector<64> data(t_block_size, 0);
            size_t pos = 1; // data[0]=0, now decode for pos > 0
            while (pos < t_block_size) {
                if (data_ptr > data_ptr_end or (data_ptr == data_ptr_end and offset >= offset_end)) {
//                    std::cout<<"entering corner case"<<std::endl;
                    while (pos < t_block_size) {
                        data[pos] = data[pos-1]+1;
//                        std::cout<<"writing entry "<<pos<<" data[pos]="<<data[pos]<<std::endl;
                        ++pos;
                    }
                    return data;
                }
                uint64_t delta = t_coder::decode(data_ptr,offset);
                if (delta == 1) {   // encoded run of ones of length >= 1
                    uint64_t len = t_coder::decode(data_ptr, offset);
                    for (size_t i=0; i<len; ++i) {
                        data[pos] = data[pos-1]+1;
                        ++pos;
                    }
                } else { // single delta
                    data[pos] = data[pos-1]+delta;
                    ++pos;
                }
            }
            return data;
        }

        static size_type estimate_size(size_type, int_vector<64>& data)
        {
            return encode(data);
        }

        static size_type
        serialize(bit_vector& bv, size_type offset, int_vector<64>& data, size_type)
        {
            return encode(data, &bv, offset);
        }

//        static size_type select_1(const bit_vector& bv, size_type offset, size_type i, size_type)
        static size_type
        select_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type)
        {
            auto offset_begin = block_start[block_id];
            auto offset_end = block_start[block_id+1];
            return decode(bv.data()+(offset_begin/64), offset_begin%64, bv.data()+(offset_end/64), offset_end%64)[i];
        }

        static size_type rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type)
        {
            auto abs_offset = block_start[block_id];
            auto data_ptr = bv.data()+(abs_offset/64);
            uint8_t offset = abs_offset%64;
            auto abs_offset_end = block_start[block_id+1];
            auto data_ptr_end = bv.data()+(abs_offset_end/64);
            uint8_t offset_end = abs_offset_end%64;

            uint64_t data_res = 0;
            size_t res = 0; // data[0]=0, now decode for pos > 0
            while (res < t_block_size and i > data_res) {
                if (data_ptr > data_ptr_end or (data_ptr == data_ptr_end and offset >= offset_end)) {
                    uint64_t len = t_block_size - res;
                    data_res += len;
                    if (i > data_res) {
                        res += len;
                    } else { // i <= data_res and i > data_res - len
                        uint64_t gap = i - (data_res - len);
                        res += gap;
                    }
                } else {
                    uint64_t delta = t_coder::decode(data_ptr,offset);
                    if (delta == 1) {   // encoded run of ones of length >= 1
                        uint64_t len = t_coder::decode(data_ptr, offset);
                        data_res += len;
                        if (i > data_res) {
                            res += len;
                        } else { // i <= data_res and i > data_res - len
                            uint64_t gap = i - (data_res - len);
                            res += gap;
                        }
                    } else { // single delta
                        data_res += delta;
                        ++res;
                    }
                }
            }
//            auto check_res = trivial_rank_1(bv, block_start, block_id, i, 0);
//            if ( check_res != res ) {
//                std::cout<<"block_id="<<block_id<<" res="<<res<<" check_res="<<check_res<<std::endl;
//            }
            return res;
        }

        static std::array<size_type,2>
        rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, const std::array<size_type,2>& ij, size_type)
        {
            auto abs_offset = block_start[block_id];
            auto data_ptr = bv.data()+(abs_offset/64);
            uint8_t offset = abs_offset%64;
            auto abs_offset_end = block_start[block_id+1];
            auto data_ptr_end = bv.data()+(abs_offset_end/64);
            uint8_t offset_end = abs_offset_end%64;

            uint64_t data_res = 0;
            std::array<size_type,2> res = {0,0}; // data[0]=0, now decode for pos > 0
            size_t k=0;
            for (; k<2; ++k) {
                while (res[k] < t_block_size and ij[k] > data_res) {
                    if (data_ptr > data_ptr_end or (data_ptr == data_ptr_end and offset >= offset_end)) {
                        uint64_t len = t_block_size - res[k];
                        data_res += len;
                        if (ij[k] > data_res) {
                            res[k] += len;
                        } else { // ij[k] <= data_res and i > data_res - len
                            uint64_t gap = ij[k] - (data_res -len);
                            res[k] += gap;
                            if (k == 0) {
                                if (ij[1] <= data_res) {
                                    res[1] = res[0] + (ij[1]-ij[0]);
                                    k = 3; break;
                                } else {
                                    res[1] = res[0] - gap + len;
                                }
                            }
                        }
                    } else {
                        uint64_t delta = t_coder::decode(data_ptr,offset);
                        if (delta == 1) {   // encoded run of ones of length >= 1
                            uint64_t len = t_coder::decode(data_ptr, offset);
                            data_res += len;
                            if (ij[k] > data_res) {
                                res[k] += len;
                            } else { // ij[k] <= data_res and i > data_res-len
                                uint64_t gap = ij[k] - (data_res-len);
                                res[k] += gap;
                                if (k == 0) {
                                    if (ij[1] <= data_res) {
                                        res[1] = res[0] + (ij[1]-ij[0]);
                                        k = 3; break;
                                    } else {
                                        res[1] = res[0] - gap + len;
                                    }
                                }
                            }
                        } else { // single delta
                            data_res += delta;
                            ++res[k];
                        }
                    }
                }
                if (k==0) {
                    res[1]= std::max(res[0], res[1]);
                }
            }
            /*
                        std::array<size_type,2> check = { rank_1(bv, block_start, block_id, ij[0], 0),
                                                          rank_1(bv, block_start, block_id, ij[1], 0)};
                        if ( res != check ){
                            std::cerr<<"res!=check"<<std::endl;
                            std::cout<<"res=["<<res[0]<<","<<res[1]<<"] != ["<<check[0]<<","<<check[1]<<"] k="<<k << std::endl;
            //                throw std::logic_error("check failed");
                            return check;
                        }
            */
            return res;
        }



        static size_type
        trivial_rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type)
//        rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type)
        {
            auto offset_begin = block_start[block_id];
            auto offset_end = block_start[block_id+1];
            auto data = decode(bv.data()+(offset_begin/64), offset_begin%64, bv.data()+(offset_end/64), offset_end%64);
            size_type res = 0;
            while (res < data.size() and i > data[res])
                ++res;
            return res;
        }

        /*
                static std::array<size_type,2>
                rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, const std::array<size_type,2>& ij, size_type)
                {
                    auto offset = block_start[block_id];
                    auto data = decode(bv.data()+(offset/64), offset%64);
                    std::array<size_type,2> res = {0,0};
                    while ( res[0] < data.size() and ij[0] > data[res[0]] )
                        ++res[0];
                    res[1] = res[0];
                    while ( res[1] < data.size() and ij[1] > data[res[1]] )
                        ++res[1];
                    return res;
                }
        */
};



template <uint64_t t_block_size>
class hyb_sd_block_full
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
        serialize(bit_vector&, size_type, int_vector<64>&, size_type u)
        {
            if (t_block_size != u) {
                std::cerr << "this should not happen!" << std::endl;
            }
            return 0;
        }

        static size_type select_1(const bit_vector&, size_type, size_type i, size_type)
        {
            return i;
        }

        static size_type rank_1(const bit_vector&, const int_vector<>&, size_type, size_type i, size_type)
        {
            return std::min(t_block_size, i);
        }
};

template <uint64_t t_block_size>
class hyb_sd_block_ef
{
    public:
        typedef bit_vector::size_type size_type;
        typedef size_type value_type;

    private:
        // TODO factor out calculation of logu and logm

    public:
        static size_type estimate_size(size_type u)
        {
            uint8_t logu = bits::hi(u) + 1;
            uint8_t logm = bits::hi(t_block_size) + 1; // TODO constexpr for hi?
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;
            size_type size_in_bits = width_low * t_block_size + (1ULL << logm) + t_block_size + 1;
            return size_in_bits;
        }
        static size_type
        serialize(bit_vector& bv, size_type offset, int_vector<64>& data, size_type u)
        {
            size_type written_bits = 0;
            uint8_t logu = bits::hi(u) + 1;
            uint8_t logm = bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;
            /* write low */
            auto data_ptr = bv.data() + (offset / 64);
            uint8_t in_word_offset = offset % 64;
            for (size_type i = 0; i < t_block_size; i++) {
                uint64_t x = data[i];
                bits::write_int_and_move(data_ptr, x, in_word_offset, width_low);
            }
            written_bits += width_low * t_block_size;

            /* write high */
            size_type last_high = 0;
            for (size_type i = 0; i < t_block_size; i++) {
                uint64_t x = data[i];
                size_type cur_high = x >> width_low;
                size_type write_val = cur_high - last_high;
                while (write_val >= 64) {
                    bits::write_int_and_move(data_ptr, 0ULL, in_word_offset, 64);
                    write_val -= 64;
                    written_bits += 64;
                }
                bits::write_int_and_move(data_ptr, 1ULL << write_val, in_word_offset, write_val + 1);
                last_high = cur_high;
                written_bits += write_val + 1;
            }
            bv[offset+written_bits] = 0;
            ++written_bits;
            return written_bits;
        }

        static size_type select_1(const bit_vector& bv, size_type offset, size_type i, size_type u)
        {
            uint8_t logu = bits::hi(u) + 1;
            uint8_t logm = bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;
            size_type hi_part_offset = offset + t_block_size * width_low;
            size_type low_part_offset = offset + i * width_low;

            auto low_part_data_ptr = bv.data() + (low_part_offset / 64);
            uint8_t low_part_in_word_offset = low_part_offset % 64;
            auto low_part = bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low);

            auto bucket = sel(bv.data(), hi_part_offset, i + 1) - hi_part_offset - i;
            return (bucket << width_low) | low_part;
        }

        static size_type rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, size_type i, size_type u)
        {
//std::cout<<">>>>>>>>rank_1("<<i<<")!!!"<<std::endl;
            auto offset = block_start[block_id];
            auto next_offset = block_start[block_id + 1];

            uint8_t logu = bits::hi(u) + 1;
            uint8_t logm = bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;

            size_type hi_part_offset = offset + t_block_size * width_low;
            size_type hi_size = next_offset - hi_part_offset;

            size_type high_val = (i >> width_low);
            size_type zeros_in_high = hi_size - t_block_size;
            if (zeros_in_high < high_val+1) {
                return t_block_size;
            }
            size_type local_sel = sel0(bv.data(), hi_part_offset, high_val + 1);

            size_type sel_high = local_sel;
            size_type rank_low = sel_high - high_val;
            if (0 == rank_low) {
                return 0;
            }

            size_type low_part_offset = offset + rank_low * width_low;
            size_type val_low = i & bits::lo_set[width_low];
            auto low_part_data_ptr = bv.data() + (low_part_offset / 64);
            uint8_t low_part_in_word_offset = low_part_offset % 64;

//std::cout<<"_sel_high="<<sel_high<<" rank_low="<<rank_low<<std::endl;
            do {
                if (!sel_high)
                    return 0;
                --rank_low;
                --sel_high;
                low_part_offset -= width_low;
                low_part_data_ptr = bv.data() + (low_part_offset / 64);
                low_part_in_word_offset = low_part_offset % 64;
//std::cout<<">>sel_high "<<bv[hi_part_offset+sel_high];
//if ( bv[hi_part_offset+sel_high] ) {
//    std::cout<<" i="<<i<<" "<< (bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low)|(high_val<<width_low))
//             <<" rank_low="<<rank_low;
//}
//std::cout<<std::endl;
            } while (bv[hi_part_offset + sel_high] and bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low) >= val_low);
            return rank_low + 1;
        }

        static std::array<size_type,2>
        rank_1(const bit_vector& bv, const int_vector<>& block_start, size_type block_id, std::array<size_type,2> ij, size_type u)
        {
            auto start_offset = block_start[block_id];
            auto next_offset = block_start[block_id + 1];

            uint8_t logu = bits::hi(u) + 1;
            uint8_t logm = bits::hi(t_block_size) + 1;
            if (logm == logu)
                logm--;
            size_type width_low = logu - logm;

            size_type hi_part_offset = start_offset + t_block_size * width_low;
            size_type hi_size = next_offset - hi_part_offset;

            std::array<size_type,2> high_val = {(ij[0] >> width_low),(ij[1] >> width_low)};

            size_type zeros_in_high = hi_size - t_block_size;
//std::cout<<"zeros_in_high="<<zeros_in_high<<std::endl;
//std::cout<<"hi_size="<<hi_size<<std::endl;
//std::cout<<"high_val[0]+1="<<high_val[0]+1<<std::endl;
//std::cout<<"high_val[1]+1="<<high_val[1]+1<<std::endl;
            if (zeros_in_high < high_val[0]+1) {    // check if there is a zero to select
                return {t_block_size, t_block_size};
            }
            std::array<size_type,2> res = {0,0};
            std::array<size_type,2> local_sel;
            local_sel[0]= sel0(bv.data(), hi_part_offset, high_val[0] + 1);
            if (high_val[0] == high_val[1]) {
                local_sel[1] = local_sel[0];
            } else { // now high_val[0] < high_val[1]
                if (zeros_in_high < high_val[1]+1) {
                    res  = {0, t_block_size}; // initialized second result
                } else {
                    if (zeros_in_high < high_val[1]+1) {    // check if there is a zero to select
                        res = {0, t_block_size};
                    } else { // there is something to select ;)
                        size_type skip = local_sel[0]+1;
//std::cout<<"skip="<<skip<<std::endl;
                        local_sel[1] = sel0(bv.data(), hi_part_offset+skip, high_val[1]-high_val[0]) + skip;
                    }
                }
            }

            bool done1 = (res[1]==t_block_size);
//            std::cout<<"done1="<<done1<<std::endl;
            for (size_t k=1-done1, s=done1; s<2; ++s, --k) {
                size_type sel_high = local_sel[k];
                size_type rank_low = sel_high - high_val[k];
//                std::cout<<"k="<<k<<" sel_high="<<sel_high<<" rank_low="<<rank_low<<std::endl;
                if (0 == rank_low) {
                    return {0,res[1]};
                }

                size_type low_part_offset = start_offset + rank_low * width_low;
                size_type val_low = ij[k] & bits::lo_set[width_low];
//std::cout<<"val_low["<<k<<"]="<<val_low<<std::endl;
                auto low_part_data_ptr = bv.data() + (low_part_offset / 64);
                uint8_t low_part_in_word_offset = low_part_offset % 64;

                do {
                    if (!sel_high) {
                        return {0, res[1]};
                    }
                    --rank_low;
                    --sel_high;
                    low_part_offset -= width_low;
                    low_part_data_ptr = bv.data() + (low_part_offset / 64);
                    low_part_in_word_offset = low_part_offset % 64;
//std::cout<<"sel_high "<<bv[hi_part_offset+sel_high];
//if ( bv[hi_part_offset+sel_high] ) {
//    std::cout<<" ij[k]="<<ij[k]<<" "<< (bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low)|(high_val[k]<<width_low))
//             <<" rank_low="<<rank_low;
//}
//std::cout<<std::endl;
                } while (bv[hi_part_offset + sel_high] and bits::read_int(low_part_data_ptr, low_part_in_word_offset, width_low) >= val_low);
                res[k] = rank_low+1;
//std::cout<<"res["<<k<<"]="<<res[k]<<std::endl;
            }
            return res;
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
        bit_vector m_bottom;
        int_vector<> m_block_start;
        int_vector<2> m_block_type;

        size_type m_size = 0;
        size_type m_num_ones = 0;


    public:
        static constexpr uint16_t block_size = t_block_size;

    private:
        enum class hyb_sd_blocktype : uint8_t {
            EF=0,
            BV=1,
            FULL=2,
            RL=3
        };

        std::pair<hyb_sd_blocktype, size_type>
        determine_block_type(size_t u, int_vector<64>& data) const
        {
            if (u == t_block_size or data[t_block_size-1] == t_block_size-1)
                return { hyb_sd_blocktype::FULL, 0 };
            std::vector<std::pair<size_type, hyb_sd_blocktype>> size_and_type;
            size_and_type.push_back({hyb_sd_block_ef<t_block_size>::estimate_size(u), hyb_sd_blocktype::EF});
            size_and_type.push_back({hyb_sd_block_bv<t_block_size>::estimate_size(u), hyb_sd_blocktype::BV});
            size_and_type.push_back({hyb_sd_block_rl<t_block_size>::estimate_size(u, data), hyb_sd_blocktype::RL});
            std::sort(size_and_type.begin(), size_and_type.end());

            /*            auto mini = std::min_element(size_and_type.begin(), size_and_type.end());
                        return { std::get<1>(*mini), std::get<0>(*mini) };
            */
            if (std::get<1>(size_and_type[0]) == hyb_sd_blocktype::RL) {
                auto rl_size = std::get<0>(size_and_type[0]);
                auto next_size = std::get<0>(size_and_type[1]);
                if (2*rl_size > next_size) {
                    std::swap(size_and_type[0], size_and_type[1]);
                }
            }
            return {std::get<1>(size_and_type[0]), std::get<0>(size_and_type[0])};
        }

        size_type compress_block(size_type i, int_vector<64>& data, size_t u)
        {
            size_type offset = m_block_start[i];
            size_type written_bits = 0;
            auto bt = determine_block_type(u, data);
            auto type = bt.first;
            auto size_in_bits = bt.second;
            m_block_type[i] = static_cast<uint8_t>(type);
            if (m_bottom.size() < offset + size_in_bits) {
                m_bottom.resize(m_bottom.size() * 2 + size_in_bits);
            }
            switch (type) {
                case hyb_sd_blocktype::BV:
                    written_bits = hyb_sd_block_bv<t_block_size>::serialize(m_bottom, offset, data, u);
                    break;
                case hyb_sd_blocktype::EF:
                    written_bits = hyb_sd_block_ef<t_block_size>::serialize(m_bottom, offset, data, u);
                    break;
                case hyb_sd_blocktype::FULL:
                    /* nothing to store */
                    break;
                case hyb_sd_blocktype::RL:
                    written_bits = hyb_sd_block_rl<t_block_size>::serialize(m_bottom, offset, data, u);
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
                        //top_lvl.push_back(i-t_block_size*top_lvl.size());
                        one_found = 0;
                    }
                }
            }
            // terminate the top level so top[i+1] - top[i] always works
//        top_lvl.push_back(last_one + num_dummy + 1);
            top_lvl.push_back(std::max(bv.size(), last_one+1) + num_dummy);

            // (2) bottom level
            m_block_start.resize(num_blocks + 1);
            m_block_type.resize(num_blocks);
            size_type value_offset = 0;
            size_type written_bits = 0;
            int_vector<64> tmp_data(t_block_size);
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
                    } else {
                        tmp_data[one_found] = std::max(m_size-value_offset, tmp_data[one_found - 1] + 1);
                        one_found++;
                    }
                    ++j;
                }
                // (2b) compress block
                size_type block_universe = top_lvl[i + 1] - top_lvl[i];
                auto wb = compress_block(i, tmp_data, block_universe);
                written_bits += wb;
            }
            m_block_start[num_blocks] = written_bits;
            m_bottom.resize(written_bits);

            // (3) encode the top level
            m_top = decltype(m_top)(top_lvl.begin(), top_lvl.end());
            m_top_sel = decltype(m_top_sel)(&m_top);
            m_top_rank = decltype(m_top_rank)(&m_top);

            // (4) bit compress pointers
            util::bit_compress(m_block_start);
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
            m_block_type.resize(num_blocks);
//std::cout<<"num_blocks="<<num_blocks<<std::endl;
            itr = begin;
            size_type value_offset = 0;
            size_type written_bits = 0;
            int_vector<64> tmp_data(t_block_size);
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
                auto wb = compress_block(i, tmp_data, block_universe);
                written_bits += wb;
            }
            m_block_start[num_blocks] = written_bits;
            m_bottom.resize(written_bits);

            // (3) encode the top level
            m_top = decltype(m_top)(top_lvl.begin(), top_lvl.end());
            m_top_sel = decltype(m_top_sel)(&m_top);
            m_top_rank = decltype(m_top_rank)(&m_top);

            // (4) bit compress pointers
            util::bit_compress(m_block_start);
//std::cout<<"bye hyb_sd"<<std::endl;
        }


        value_type operator[](size_type i)const
        {
            auto ranks = rank_1({{i+1,i}});
            return ranks[1]-ranks[0];
//            return rank_1(i+1) - rank_1(i);
        }

        //! Accessing the i-th element of the original bit_vector
        size_type select_1(size_type i) const
        {
            i = i - 1;
            auto block_id = i / t_block_size;
            auto in_block_offset = i % t_block_size;
            auto top_value = m_top_sel(block_id + 1);
            size_type res = top_value;

            if (in_block_offset == 0)
                return res;

            auto u = m_top_sel(block_id + 2) - top_value;
            /*
            auto bt = determine_block_type(u);
            auto block_type = bt.first;
            */

            auto block_type = static_cast<hyb_sd_blocktype>(m_block_type[block_id]);
            size_type block_offset = m_block_start[block_id];

            switch (block_type) {
                case hyb_sd_blocktype::BV:
                    res += hyb_sd_block_bv<t_block_size>::select_1(m_bottom, block_offset, in_block_offset, u);
                    break;
                case hyb_sd_blocktype::EF:
                    res += hyb_sd_block_ef<t_block_size>::select_1(m_bottom, block_offset, in_block_offset, u);
                    break;
                case hyb_sd_blocktype::FULL:
                    res += hyb_sd_block_full<t_block_size>::select_1(m_bottom, block_offset, in_block_offset, u);
                    break;
                case hyb_sd_blocktype::RL:
                    res += hyb_sd_block_rl<t_block_size>::select_1(m_bottom, m_block_start, block_id, in_block_offset, u);
                    break;
            }
            return res;
        }

        size_type rank_1(size_type i) const
        {
//            bool debug = false;
//            if ( i==2075 or i==2076) {
//                debug = true;
//            }
//std::cout<<"!!! rank_1("<<i<<")"<<std::endl;
            if (i > m_size or m_num_ones == 0) {
//std::cout<<"!!! i > m_size "<<i<<" > "<<m_size<<std::endl;
                return m_num_ones;
            }
            auto block_id = m_top_rank(i);
            if (block_id == 0) {
//std::cout<<"!!! block_id=0"<<std::endl;
                return 0;
            }
            block_id -= 1;
            size_type res = block_id * t_block_size;
            auto top_value = m_top_sel(block_id + 1);
            size_type in_block_i = i;
            in_block_i -= top_value;
            if (in_block_i == 0) {
//std::cout<<"!!! in_block_i=0"<<std::endl;
                return res;
            }
// TODO: can we return res+in_block_i if top_value-i==in_block_i ???

            auto block_type = static_cast<hyb_sd_blocktype>(m_block_type[block_id]);
//if (debug) {
//    std::cout<<"!!! i="<<i<<" block_id="<<block_id<<" res="<<res<<std::endl;
//    std::cout<<"!!! hyb_sd_blocktype = "<< (int)block_type <<std::endl;
//}
            if (block_type == hyb_sd_blocktype::FULL) {
//if (debug) std::cout<<"!!! hyb_sd_blocktype::FULL"<<std::endl;
                return res + hyb_sd_block_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, 0);
            }

            auto u = m_top_sel(block_id + 2) - top_value;
//            auto bt = determine_block_type(u);
//            auto block_type = bt.first;
            size_type block_offset = m_block_start[block_id];

            switch (block_type) {
                case hyb_sd_blocktype::BV:
//             std::cout << "!!!BV" << std::endl;
                    res += hyb_sd_block_bv<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
                case hyb_sd_blocktype::EF:
//if (debug){
//    std::cout << "!!!single EF in_block_i="<<in_block_i << std::endl;
//}
                    res += hyb_sd_block_ef<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
                case hyb_sd_blocktype::FULL:
//             std::cout << "!!!FULL" << std::endl;
                    res += hyb_sd_block_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
                case hyb_sd_blocktype::RL:
//             std::cout << "!!!RL" << std::endl;
                    res += hyb_sd_block_rl<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    break;
            }
//if (debug) {
//    std::cout<<"!!!!  res="<<res<<std::endl;
//}
            return res;
        }

        std::array<size_type,2>
        rank_1(std::array<size_type,2> ij) const
        {
            if (ij[0] > ij[1]) {
                return {rank_1(ij[0]),rank_1(ij[1])};
            }
            // no we know ij[0] <= ij[1]
            if (ij[0] > m_size or m_num_ones == 0) {
                return {m_num_ones, m_num_ones};
            }
            if (ij[1] > m_size or m_num_ones == 0) {
                return {rank_1(ij[0]), m_num_ones};
            }
            auto block_id = m_top_rank(ij[0]);
            if (block_id == 0) {
                size_type first_element = m_top_sel(1);
                if (ij[1] <= first_element) {
                    return {0,0};
                }
                return {0, rank_1(ij[1])}; // TODO: can still be optimized
            }
            block_id -= 1;
            size_type r = block_id * t_block_size;
            auto top_value = m_top_sel(block_id + 1);
            size_type in_block_i = ij[0];
            in_block_i -= top_value;
            size_type in_block_j = ij[1];
            in_block_j -= top_value;

            if (in_block_i == 0) {
                if (ij[0]==ij[1]) {
                    return {r,r};
                }
                return {r, rank_1(ij[1])}; // TODO: can still be optimized
            }

            auto block_type = static_cast<hyb_sd_blocktype>(m_block_type[block_id]);
            if (block_type == hyb_sd_blocktype::FULL and in_block_j < t_block_size) {
                return {r+hyb_sd_block_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, 0),
                        r+hyb_sd_block_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_j, 0)
                       };
            }

            auto u = m_top_sel(block_id + 2) - top_value;
//            auto bt = determine_block_type(u);
//            auto block_type = bt.first;
            size_type block_offset = m_block_start[block_id];
            std::array<size_type, 2> res {r,r};

            switch (block_type) {
                case hyb_sd_blocktype::BV:
//                    std::cout<<"double rank_1 for BV"<<std::endl;
                    if (in_block_j >= u) {
                        res[0] += hyb_sd_block_bv<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    } else {
                        auto in_block_rank = hyb_sd_block_bv<t_block_size>::rank_1(m_bottom, m_block_start, block_id, {in_block_i, in_block_j}, u);
                        res[0] += in_block_rank[0];
                        res[1] += in_block_rank[1];
                    }
                    break;
                case hyb_sd_blocktype::EF:
//                    std::cout<<"double rank_1 for EF"<<std::endl;
                    if (in_block_j >= u) {
                        res[0] += hyb_sd_block_ef<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    } else {
//                        std::cout<<"call double"<<std::endl;
                        auto in_block_rank = hyb_sd_block_ef<t_block_size>::rank_1(m_bottom, m_block_start, block_id, {in_block_i, in_block_j}, u);
                        res[0] += in_block_rank[0];
                        res[1] += in_block_rank[1];
                    }
                    break;
                case hyb_sd_blocktype::FULL:
//                    std::cout<<"double rank_1 for FULL"<<std::endl;
                    res[0] += hyb_sd_block_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    if (in_block_j < u) {
                        res[1] += hyb_sd_block_full<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_j, u);
                    }
                    break;
                case hyb_sd_blocktype::RL:
//                    std::cout<<"double rank_1 for RL"<<std::endl;
                    if (in_block_j >= u) {
                        res[0] += hyb_sd_block_rl<t_block_size>::rank_1(m_bottom, m_block_start, block_id, in_block_i, u);
                    } else {
                        auto in_block_rank = hyb_sd_block_rl<t_block_size>::rank_1(m_bottom, m_block_start, block_id, {in_block_i, in_block_j}, u);
                        res[0] += in_block_rank[0];
                        res[1] += in_block_rank[1];
                    }
                    break;
            }
            if (in_block_j >= u) {
                res[1] = rank_1(ij[1]);
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
            //written_bytes += m_bottom.serialize(out, nullptr, "bottom");
            auto bottom_bytes = m_bottom.serialize(out, nullptr, "bottom");
            {
                structure_tree_node* bottom_child = structure_tree::add_child(child, "bottom", util::class_name(m_bottom));
                std::array<size_type,4> written_bits = {0,0,0,0};
                for (size_t i=1; i<m_block_start.size(); ++i) {
                    written_bits[m_block_type[i-1]] += m_block_start[i]-m_block_start[i-1];
                }
                std::vector<std::string> names = {"EF","BV","FULL","RL"};
                for (size_t i=0; i<written_bits.size(); ++i) {
                    structure_tree_node* block_child = structure_tree::add_child(bottom_child, names[i], util::class_name(m_bottom));
                    structure_tree::add_size(block_child, (written_bits[i]+7)/8);
                }
                structure_tree::add_size(bottom_child, bottom_bytes);
            }
            written_bytes += bottom_bytes;
            written_bytes += m_block_start.serialize(out, child, "block_start");
            written_bytes += m_block_type.serialize(out, child, "block_type");
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
            m_block_type.load(in);
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
            m_block_type.swap(v.m_block_type);
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

        template<typename t_pos>
        t_pos rank(t_pos i) const
        {
            return m_v->rank_1(i);
        }

        template<typename t_pos>
        t_pos operator()(t_pos i) const
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
