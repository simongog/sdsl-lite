#pragma once

#include <iterator>

#include "bit_streams.hpp"
#include "bit_coders.hpp"
#include "bit_magic.hpp"

#include "list_basics.hpp"

template<uint16_t t_skip,bool t_sorted,bool t_compact>
class ef_skip_iterator : public std::iterator<std::random_access_iterator_tag,uint64_t,std::ptrdiff_t>
{
    public:
        using size_type = sdsl::int_vector<>::size_type;
        uint64_t SKIP_OFFSET_THRESHOLD = 16;
    private:
        uint64_t m_size = 0;
        uint64_t m_universe = 0;
        uint8_t m_width_low = 0;
        uint64_t m_low_offset = 0;
        uint64_t m_high_offset = 0;
        uint64_t m_start_high = 0;
        uint64_t m_skip_width = 0;
        uint64_t m_skip_start_offset = 0;
        const uint64_t* m_data;
        bool m_end = false;
    private:
        mutable value_type m_cur_elem = 0;
        mutable size_type m_last_accessed_offset = std::numeric_limits<uint64_t>::max();
        size_type m_cur_high_offset = 0;
        size_type m_cur_offset = 0;
        mutable value_type m_prev_elem = 0;
    public:
        ef_skip_iterator(const bit_istream& is,size_t start_offset,bool end) : m_end(end)
        {
            static_assert(t_compact == false,"Can't use this constructor with compact ef representation.");
            m_data = is.data();
            is.seek(start_offset);
            m_size = is.decode<coder::elias_gamma>();
            m_universe = is.decode<coder::elias_gamma>();
            uint8_t logm = sdsl::bits::hi(m_size)+1;
            uint8_t logu = sdsl::bits::hi(m_universe)+1;
            if (logu < logm) {
                m_width_low = 1;
            } else {
                if (logm == logu) logm--;
                m_width_low = logu - logm;
            }

            if (t_skip) {
                size_type num_zeros = (m_universe >> m_width_low);
                size_t num_skips = num_zeros/t_skip  + 1;
                m_skip_width = is.decode<coder::elias_gamma>();
                m_skip_start_offset = is.tellg();
                is.skip(num_skips*m_skip_width);
            }

            m_low_offset = is.tellg();
            m_high_offset = m_low_offset + m_size*m_width_low;
            m_start_high = m_low_offset + m_size*m_width_low;
            m_cur_offset = end ? m_size : 0;
            m_cur_high_offset = 0;
            if (m_cur_offset == 0) {
                if (high(0) != 1) {
                    m_cur_high_offset = sdsl::bits::next(m_data,m_high_offset+1) - m_high_offset;
                }
            }
            size_t cur_bucket = m_cur_high_offset - m_cur_offset;
            m_cur_elem = (cur_bucket << m_width_low) | low(m_cur_offset);
            m_last_accessed_offset = m_cur_offset;
        }
        ef_skip_iterator(const bit_istream& is,size_t start_offset,bool end,size_type size,size_type universe)
            : m_size(size) , m_universe(universe) , m_end(end)
        {
            uint8_t logm = sdsl::bits::hi(m_size)+1;
            uint8_t logu = sdsl::bits::hi(m_universe)+1;
            if (logu < logm) {
                m_width_low = 1;
            } else {
                if (logm == logu) logm--;
                m_width_low = logu - logm;
            }
            m_data = is.data();
            is.seek(start_offset);

            if (t_skip) {
                size_type num_zeros = (m_universe >> m_width_low);
                size_t num_skips = num_zeros/t_skip  + 1;
                m_skip_width = is.decode<coder::elias_gamma>();
                m_skip_start_offset = is.tellg();
                is.skip(num_skips*m_skip_width);
            }

            m_low_offset = is.tellg();
            m_high_offset = m_low_offset + m_size*m_width_low;
            m_start_high = m_low_offset + m_size*m_width_low;
            m_cur_offset = end ? m_size : 0;
            m_cur_high_offset = 0;
            if (m_cur_offset == 0) {
                if (high(0) != 1) {
                    m_cur_high_offset = sdsl::bits::next(m_data,m_high_offset+1) - m_high_offset;
                }
            }
            size_t cur_bucket = m_cur_high_offset - m_cur_offset;
            m_cur_elem = (cur_bucket << m_width_low) | low(m_cur_offset);
            m_last_accessed_offset = m_cur_offset;
        }
        ef_skip_iterator() = default;
        ef_skip_iterator(const ef_skip_iterator& pi) = default;
        ef_skip_iterator(ef_skip_iterator&& pi) = default;
        ef_skip_iterator& operator=(const ef_skip_iterator& pi) = default;
        ef_skip_iterator& operator=(ef_skip_iterator&& pi) = default;
    public:
        size_type offset() const
        {
            return m_cur_offset;
        }
        size_type size() const
        {
            return m_size;
        }
        size_type remaining() const
        {
            return m_size - m_cur_offset;
        }
        uint64_t operator*() const
        {
            if (!t_sorted) {
                if (m_last_accessed_offset != m_cur_offset) {
                    if (m_cur_offset == 0 || m_last_accessed_offset+1 == m_cur_offset) {
                        m_prev_elem = m_cur_elem;
                    } else {
                        size_type poffset = m_high_offset + m_cur_high_offset;
                        auto prev_high_offset = sdsl::bits::prev(m_data,poffset-1) - m_high_offset;
                        size_t prev_bucket = prev_high_offset - (m_cur_offset - 1);
                        m_prev_elem = (prev_bucket << m_width_low) | low(m_cur_offset-1);
                        size_t cur_bucket = m_cur_high_offset - m_cur_offset;
                        m_cur_elem = (cur_bucket << m_width_low) | low(m_cur_offset);
                    }
                    size_t cur_bucket = m_cur_high_offset - m_cur_offset;
                    m_cur_elem = (cur_bucket << m_width_low) | low(m_cur_offset);
                    m_last_accessed_offset = m_cur_offset;
                }
                return m_cur_elem - m_prev_elem;
            } else {
                if (m_last_accessed_offset != m_cur_offset) {
                    size_t cur_bucket = m_cur_high_offset - m_cur_offset;
                    m_cur_elem = (cur_bucket << m_width_low) | low(m_cur_offset);
                    m_last_accessed_offset = m_cur_offset;
                }
                return m_cur_elem;
            }
        }
        bool operator ==(const ef_skip_iterator& b) const
        {
            return (b.m_end && m_cur_offset > b.m_cur_offset) || m_cur_offset == b.m_cur_offset;
        }
        bool operator !=(const ef_skip_iterator& b) const
        {
            return m_cur_offset != b.m_cur_offset;
        }
        ef_skip_iterator operator++(int)
        {
            ef_skip_iterator tmp(*this);
            size_type offset = m_high_offset + m_cur_high_offset;
            m_cur_high_offset = sdsl::bits::next(m_data,offset+1) - m_high_offset;
            m_cur_offset++;
            return tmp;
        }
        ef_skip_iterator& operator++()
        {
            size_type offset = m_high_offset + m_cur_high_offset;
            m_cur_high_offset = sdsl::bits::next(m_data,offset+1) - m_high_offset;
            m_cur_offset++;
            return *this;
        }
        ef_skip_iterator& operator+=(size_type i)
        {
            if (i == 0) return *this;
            size_type offset = m_high_offset + m_cur_high_offset;
            m_cur_high_offset = bit_magic::next_Xth_one(m_data,offset,i) - m_high_offset;
            m_cur_offset += i;
            return *this;
        }
        ef_skip_iterator operator+(size_type i)
        {
            ef_skip_iterator tmp(*this);
            tmp += i;
            return tmp;
        }
        ef_skip_iterator& operator--()
        {
            size_type offset = m_high_offset + m_cur_high_offset;
            m_cur_high_offset = sdsl::bits::prev(m_data,offset-1) - m_high_offset;
            m_cur_offset--;
            return *this;
        }
        template<class t_itr>
        auto operator-(const t_itr& b) const -> difference_type
        {
            return (difference_type)offset() - (difference_type)b.offset();
        }
        bool skip(uint64_t pos)
        {
            static_assert(t_sorted == true,"skipping only works in sorted lists.");
            if (m_cur_elem == pos) return true;
            if (m_cur_elem > pos) return false;
            if (m_universe < pos) {
                m_cur_offset = m_size;
                return false;
            }
            uint64_t high_bucket = pos >> m_width_low;
            uint64_t cur_bucket = m_cur_high_offset - m_cur_offset;
            if (cur_bucket > high_bucket) return false;
            if (high_bucket != cur_bucket) {
                // skip to the one after high_bucket 0's
                auto skip_zeros = high_bucket - cur_bucket;
                size_t zero_pos = 0;
                if (skip_zeros > SKIP_OFFSET_THRESHOLD && high_bucket > t_skip) {
                    auto high_skip = skip_offset(high_bucket);
                    auto skip_fast = high_bucket % t_skip;
                    if (skip_fast) zero_pos = bit_magic::next_Xth_zero(m_data,high_skip,skip_fast) - m_high_offset;
                    else zero_pos = high_skip - m_high_offset;
                } else {
                    zero_pos = bit_magic::next_Xth_zero(m_data,m_high_offset+m_cur_high_offset,skip_zeros) - m_high_offset;
                }
                m_cur_offset = m_cur_offset + (zero_pos - m_cur_high_offset) - skip_zeros + 1;
                m_cur_high_offset = sdsl::bits::next(m_data,m_high_offset+zero_pos) - m_high_offset;
                if (m_cur_high_offset != zero_pos+1) {
                    // if the bucket we hit doesn't exactly match we don't have to search inside the bucket
                    return false;
                }
                if (m_cur_offset > m_size) {
                    m_cur_offset = m_size;
                    return false;
                }
            }
            // scan the lower part of the current bucket
            auto next_zero_pos = bit_magic::next0(m_data,m_high_offset+m_cur_high_offset) - m_high_offset;
            auto remaining_bucket_size = next_zero_pos - m_cur_high_offset;
            size_t i = 0;
            const uint64_t low_num = pos&((1ULL<<m_width_low)-1);
            for (i=m_cur_offset; i<m_cur_offset+remaining_bucket_size; i++) {
                auto cur_low = low(i);
                if (cur_low >= low_num) {
                    // found pos or found a larger within the bucket!
                    auto ones_skipped = i-m_cur_offset;
                    m_cur_high_offset = m_cur_high_offset + ones_skipped;
                    m_cur_offset += ones_skipped;
                    if (m_cur_offset >= m_size) {
                        m_cur_offset = m_size;
                        return false;
                    }
                    m_cur_elem = (high_bucket << m_width_low) | cur_low;
                    m_last_accessed_offset = m_cur_offset;
                    if (m_cur_elem == pos) return true;
                    return false;
                }
            }
            // nothing found in the bucket we were interested in. move on to the next larger bucket
            m_cur_high_offset = sdsl::bits::next(m_data,m_high_offset+next_zero_pos) - m_high_offset;
            m_cur_offset += remaining_bucket_size;
            if (m_cur_offset >= m_size) {
                m_cur_offset = m_size;
                return false;
            }
            return false;
        }
    private:
        inline value_type low(size_type i) const
        {
            const auto off = m_low_offset + i*m_width_low;
            const auto data_ptr = m_data + (off>>6);
            const auto in_word_offset = off&0x3F;
            return sdsl::bits::read_int(data_ptr,in_word_offset,m_width_low);
        }
        inline value_type high(size_type i) const
        {
            const auto off = m_high_offset + i;
            const auto data_ptr = m_data + (off>>6);
            const auto in_word_offset = off&0x3F;
            return sdsl::bits::read_int(data_ptr,in_word_offset,1);
        }
        inline value_type skip_offset(size_type i) const
        {
            const auto skip_pos = i/t_skip;
            const auto off = m_skip_start_offset + (skip_pos*m_skip_width);
            const auto data_ptr = m_data + (off>>6);
            const auto in_word_offset = off&0x3F;
            return m_high_offset + sdsl::bits::read_int(data_ptr,in_word_offset,m_skip_width);
        }
};

template<uint16_t t_skip = 64,bool t_sorted = true,bool t_compact = false>
struct eliasfano_skip_list {
    using size_type = sdsl::int_vector<>::size_type;
    using iterator_type = ef_skip_iterator<t_skip,t_sorted,t_compact>;
    using list_type = list_dummy<iterator_type>;

    template<class t_itr>
    static size_type create(bit_ostream& os,t_itr begin,t_itr end,size_type m = 0,size_type u = 0)
    {
        size_type data_offset = os.tellp();

        // compute properties
        if (m == 0) m = std::distance(begin,end);
        if (u == 0) {
            if (!t_sorted) {
                u = std::accumulate(begin,end, 0LL)+1;
            } else {
                u = *(end-1)+1;
            }
        }
        uint8_t logm = sdsl::bits::hi(m)+1;
        uint8_t logu = sdsl::bits::hi(u)+1;
        uint8_t width_low;
        if (logu < logm) {
            width_low = 1;
        } else {
            if (logm == logu) logm--;
            width_low = logu - logm;
        }

        if (!t_compact) {
            // write size
            os.encode<coder::elias_gamma>(m);
            // write universe size
            os.encode<coder::elias_gamma>(u);
        }

        size_type num_zeros = (u >> width_low);
        size_t skip_start_offset = 0;
        if (t_skip != 0) {
            auto num_skips = num_zeros/t_skip + 1;
            size_t skip_ptr_width = sdsl::bits::hi(num_zeros+m+1)+1;
            os.encode<coder::elias_gamma>(skip_ptr_width);
            os.expand_if_needed(num_skips*skip_ptr_width);
            skip_start_offset = os.tellp();
            os.skip(num_skips*skip_ptr_width);
        }

        // write low
        os.expand_if_needed(m*width_low);
        uint64_t last = 0;
        auto itr = begin;
        while (itr != end) {
            auto position = *itr;
            if (!t_sorted) position += last;
            last = position;
            os.put_int_no_size_check(position,width_low);
            ++itr;
        }

        // write high
        os.expand_if_needed(num_zeros+m+2);
        itr = begin;
        size_type last_high=0;
        last = 0;
        size_type start_high = os.tellp();
        while (itr != end) {
            auto position = *itr;
            if (!t_sorted) position += last;
            last = position;
            size_type cur_high = position >> width_low;
            size_type cur_num = cur_high-last_high;
            os.put_unary_no_size_check(cur_num);
            last_high = cur_high;
            ++itr;
        }
        // terminator bit so no checks for end of bitstream
        // are needed inside the iterators
        os.put(1);
        // write the skips
        if (t_skip) {
            bit_istream bis(os.bitvector(),start_high);
            size_type stop_high = os.tellp();
            size_t n_high = stop_high - start_high;
            size_t skip_ptr_width = sdsl::bits::hi(num_zeros+m+1)+1;
            auto last_pos = os.tellp();
            os.seek(skip_start_offset);
            os.put_int_no_size_check(0,skip_ptr_width);
            size_type zeros_seen = 0;
            for (size_t i=0; i<n_high; i++) {
                if (bis.get() == 0) {
                    zeros_seen++;
                    if (zeros_seen % t_skip == 0) {
                        os.put_int_no_size_check(i,skip_ptr_width);
                    }
                }
            }
            // seek back to the end of the list
            os.seek(last_pos);
        }


        return data_offset;
    }

    static list_dummy<iterator_type> materialize(const bit_istream& is,size_t start_offset)
    {
        return list_dummy<iterator_type>(iterator_type(is,start_offset,false),iterator_type(is,start_offset,true));
    }

    static list_dummy<iterator_type> materialize(const bit_istream& is,size_t start_offset,size_type m,size_type u)
    {
        return list_dummy<iterator_type>(iterator_type(is,start_offset,false,m,u),iterator_type(is,start_offset,true,m,u));
    }

    static size_type estimate_size(size_type m,size_type u)
    {
        uint8_t logm = sdsl::bits::hi(m)+1;
        uint8_t logu = sdsl::bits::hi(u)+1;
        if (logm == logu) logm--;
        uint8_t width_low = logu - logm;
        size_type num_zeros = (u >> width_low);
        size_t skip_size = 0;
        if (t_skip) {
            auto num_skips = m/t_skip;
            size_t skip_ptr_width = sdsl::bits::hi(num_zeros+m+1)+1;
            skip_size = num_skips*skip_ptr_width;
        }
        return width_low*m + (num_zeros+m+1) + skip_size;
    }
};
