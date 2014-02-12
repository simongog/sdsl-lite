/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file wt_gmr.hpp
    \brief wt_gmr.hpp contains a specialized class to support select, rank
			and access on inputs over a large alphabet.
    \author Alexander Diehm, Timo Beller, Simon Gog
*/
#ifndef INCLUDED_SDSL_WT_GMR
#define INCLUDED_SDSL_WT_GMR

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<class t_bitvector = bit_vector,
         class t_select = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
class wt_gmr_1
{
    private:

        t_bitvector m_bv_blocks;
        int_vector<> e;					//rename?
        t_select m_bv_blocks_select1;
        t_select_zero m_bv_blocks_select0;
        uint64_t m_size; // input length
        uint64_t m_max_symbol = 0; // maximum character + 1
        uint64_t m_blocks; // blocks per character
        uint64_t m_sigma =0;

    public:
        typedef int_vector<>::size_type size_type;
        typedef int_vector<>::value_type value_type;
        typedef wt_tag index_category;
        typedef int_alphabet_tag alphabet_category;
        enum {lex_ordered=0};

        const size_type&       sigma = m_sigma; // Todo

        wt_gmr_1() {}

        template<uint8_t int_width>
        wt_gmr_1(int_vector_buffer<int_width>& input, size_type size) : m_size(size) {
            // Determin max. symbol
            for (uint64_t i=0; i<m_size; ++i) {
                if (m_max_symbol < input[i]) m_max_symbol = input[i];
            }
            ++m_max_symbol;

            // Create and fill m_bv_blocks
            m_blocks = (m_size+m_max_symbol-1)/m_max_symbol;
            bit_vector b(m_size+m_max_symbol*m_blocks+1,0);
            int_vector<> symbols(m_max_symbol,0,bits::hi(m_size)+1);
            {
                int_vector<> tmp(m_max_symbol*m_blocks,0,bits::hi(m_max_symbol)+1);

                for (uint64_t i=0, offset=0, j=0; i<m_size; ++i, ++j) {
                    if (j==m_max_symbol) {
                        ++offset;
                        j = 0;
                    }
                    ++tmp[offset+input[i]*m_blocks];
                }

                for (uint64_t i=0; i<symbols.size(); ++i) {
                    for (uint64_t j=m_blocks*i; j<(i+1)*m_blocks; ++j) {
                        symbols[i] += tmp[j];
                    }
                }

                for (uint64_t i=0,l=1; i<tmp.size(); ++i,++l) {
                    for (uint64_t j=0; j<tmp[i]; ++j) b[l++]=1;
                }

                //calc m_sigma
                bool write = true;
                uint64_t blocks=0;
                for (uint64_t i=1; i<b.size(); ++i) {
                    if (blocks==m_blocks) {
                        blocks = 0;
                        write = true;
                    }
                    if (b[i]) {
                        if (write) {
                            ++m_sigma;
                            write = false;
                        }
                    } else ++blocks;
                }

                m_bv_blocks = t_bitvector(std::move(b));
            }
            util::init_support(m_bv_blocks_select0, &m_bv_blocks);
            util::init_support(m_bv_blocks_select1, &m_bv_blocks);

            // Create and fill e
            e = int_vector<>(m_size,0,bits::hi(m_max_symbol)+1);
            for (uint64_t i=0, tmp=0, sum=0; i<m_max_symbol; ++i) {
                tmp = symbols[i];
                symbols[i] = sum;
                sum += tmp;
            }
            for (uint64_t i=0; i<m_size;) {
                for (uint64_t j=0; j<m_max_symbol and i<m_size; ++i, ++j) {
                    e[symbols[input[i]]++] = j;
                }
            }
        }

        //! Swap operator
        void swap(wt_gmr_1& fs) {
            if (this != &fs) {
                m_bv_blocks.swap(fs.m_bv_blocks);
                e.swap(fs.e);
                util::swap_support(m_bv_blocks_select0, fs.m_bv_blocks_select0, &m_bv_blocks, &(fs.m_bv_blocks));
                util::swap_support(m_bv_blocks_select1, fs.m_bv_blocks_select1, &m_bv_blocks, &(fs.m_bv_blocks));
                std::swap(m_size, fs.m_size);
                std::swap(m_max_symbol, fs.m_max_symbol);
                std::swap(m_blocks, fs.m_blocks);
                std::swap(m_sigma, fs.m_sigma);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        value_type operator[](size_type i)const {
            assert(i<m_size);
            size_type block = i/m_max_symbol+1, val = i%m_max_symbol, search_begin, search_end, j;
            while (true) {
		j = m_bv_blocks_select0(block)+1;
		search_begin = j-block;
                if (m_bv_blocks[j]) {
                    search_end = m_bv_blocks_select0(block+1)-(block);
                    if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                        while (search_begin < search_end and e[search_begin] <= val) {
                            if (e[search_begin]==val) {
                                return (block-1)/m_blocks;
                            }
                            ++search_begin;
                        }
                    } else {
                        if (binary_search(e.begin()+search_begin, e.begin()+search_end, val)) {
                            return (block-1)/m_blocks;
                        }
                    }
                }
                block+=m_blocks;
            }
        }

        std::pair<size_type, value_type> inverse_select(size_type i)const {
            assert(i<m_size);
            size_type block = i/m_max_symbol+1, val = i%m_max_symbol, offset = 0, search_begin, search_end, j;
            while (true) {
                j = m_bv_blocks_select0(block)+1;
                search_begin = j-block;
                if (m_bv_blocks[j]) {
                    search_end = m_bv_blocks_select0(block+1)-(block);
                    offset = 0;
                    if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                        while (search_begin < search_end and e[search_begin] <= val) {
                            if (e[search_begin]==val) {
                                value_type c = (block-1)/m_blocks;

                                size_type ones_before_cblock = m_bv_blocks_select0(c*m_blocks+1)-(c*m_blocks);

                                size_type r = search_begin-ones_before_cblock;
                                return std::make_pair(r,c);
                            }
                            ++search_begin;
                        }
                    } else {
                        offset = lower_bound(e.begin()+search_begin, e.begin()+search_end, val)-e.begin();
                        if (offset<search_end) {
                            if (e[offset]==val) {
                                value_type c = (block-1)/m_blocks;

                                size_type ones_before_cblock = m_bv_blocks_select0(c*m_blocks+1)-(c*m_blocks);

                                size_type r = offset-ones_before_cblock;
                                return std::make_pair(r,c);
                            }
                        }
                    }
                }
                block+=m_blocks;
            }
        }

        size_type select(size_type i, value_type c)const {
            size_type k = m_bv_blocks_select0(c*m_blocks+1)-(c*m_blocks)+i;
            return (m_bv_blocks_select1(k)-k)*m_max_symbol+e[k-1]-c*m_blocks*m_max_symbol;
        }

        size_type rank(size_type i, value_type c)const {
            if (c>m_max_symbol-1) return 0;
            if (i<=0) return 0;

            size_type offset=0;
            size_type ones_before_cblock = m_bv_blocks_select0(c*m_blocks+1)-c*m_blocks;
            size_type search_begin = m_bv_blocks_select0(c*m_blocks+(i-1)/m_max_symbol+1)-(c*m_blocks+(i-1)/m_max_symbol+1)+1;
            size_type search_end = m_bv_blocks_select0(c*m_blocks+(i-1)/m_max_symbol+2)-(c*m_blocks+(i-1)/m_max_symbol+1);

            size_type val = (i-1)%m_max_symbol;
            if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                while (search_begin < search_end and e[search_begin] <= val) {
                    ++search_begin;
                }
            } else {
                offset = lower_bound(e.begin()+search_begin, e.begin()+search_end, val+1)-e.begin()-search_begin;
            }
            return search_begin+offset-ones_before_cblock;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "m_size");
            written_bytes += write_member(m_max_symbol, out, child, "m_max_symbol");
            written_bytes += write_member(m_blocks, out, child, "m_blocks");
            written_bytes += write_member(m_sigma, out, child, "m_sigma");
            written_bytes += m_bv_blocks.serialize(out, child, "m_bv_blocks");
            written_bytes += e.serialize(out, child, "e");
            written_bytes += m_bv_blocks_select0.serialize(out, child, "m_bv_blocks_select0");
            written_bytes += m_bv_blocks_select1.serialize(out, child, "m_bv_blocks_select1");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_max_symbol, in);
            read_member(m_blocks, in);
            read_member(m_sigma, in);
            m_bv_blocks.load(in);
            e.load(in);
            m_bv_blocks_select0.load(in, &m_bv_blocks);
            m_bv_blocks_select1.load(in, &m_bv_blocks);
        }
};

template<class t_bitvector = bit_vector,
         class t_select = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_rank = typename t_bitvector::rank_1_type>
class wt_gmr_2
{
    private:

        t_bitvector m_bv_blocks; // 0 indicates end of block. Corresponds to B in the paper.
        t_bitvector m_bv_chunks; // 0 indicates end of symbol in chunk. Corresponds to X in the paper.
        t_bitvector m_bv_shortcut; // 1 indicates a shortcut

        int_vector<> m_shortcut; // Contains the shortcuts. Corresponds to S in the paper.
        int_vector<> m_perm; // Contains permutation of each chunk. Corresponds to \f$ \pi \f$ in the paper.

        t_rank m_bv_shortcut_rank;
        t_select m_bv_blocks_select1, m_bv_chunks_select1;
        t_select_zero m_bv_blocks_select0, m_bv_chunks_select0;

        uint64_t m_size; // input length
        uint64_t m_t = 32; // shortcut value todo: set via template or construct
        uint64_t m_max_symbol = 0; // maximum character + 1
        uint64_t m_chunks; // number of chunks

        uint64_t m_sigma = 0;

    public:
        typedef int_vector<>::size_type size_type;
        typedef int_vector<>::value_type value_type;
        typedef wt_tag index_category;
        typedef int_alphabet_tag alphabet_category;
        enum {lex_ordered=0};

        const size_type&       sigma = m_sigma; // Todo

        wt_gmr_2() {}

        template<uint8_t int_width>
        wt_gmr_2(int_vector_buffer<int_width>& input, size_type size) : m_size(size) {
            // Determin max. symbol
            for (uint64_t i=0; i<m_size; ++i) {
                if (m_max_symbol < input[i]) m_max_symbol = input[i];
            }
            ++m_max_symbol;

            m_chunks = (m_size+m_max_symbol-1)/m_max_symbol;


            m_perm = int_vector<>(m_size,0,bits::hi(m_max_symbol-1)+1);
            {
                uint64_t x_pos = 0;
                bit_vector x(m_max_symbol*m_chunks+m_size+1,0);

                //fill m_perm and m_bv_chunks for every chunk
                for (uint64_t i=0; i<m_chunks; ++i) {
                    int_vector<> symbols(m_max_symbol,0,bits::hi(m_max_symbol-1)+2);

                    //calc symbols
                    for (uint64_t j=i*m_max_symbol; j<(i+1)*m_max_symbol and j<m_size; ++j) {
                        ++symbols[input[j]];
                    }
                    //calc m_bv_chunks
                    for (uint64_t j=0; j<m_max_symbol; ++j,++x_pos) {
                        for (uint64_t k=0; k<symbols[j]; ++k) x[++x_pos]=1;
                    }

                    //calc symbols prefix sum
                    for (uint64_t j=0, tmp=0, sum=0; j<m_max_symbol; ++j) {
                        tmp = symbols[j];
                        symbols[j] = sum;
                        sum += tmp;
                    }
                    //calc m_perm
                    for (uint64_t j=i*m_max_symbol, k=0; j<(i+1)*m_max_symbol and j<m_size; ++j,++k) {
                        m_perm[i*m_max_symbol+(symbols[input[j]]++)] = k;
                    }
                }
                m_bv_chunks = t_bitvector(std::move(x));
                util::init_support(m_bv_chunks_select1, &m_bv_chunks);
                util::init_support(m_bv_chunks_select0, &m_bv_chunks);

            }

            //calc m_bv_blocks
            {
                bit_vector b(m_size+m_max_symbol*m_chunks+1,0);
                int_vector<> tmp(m_max_symbol*m_chunks,0,bits::hi(m_max_symbol-1)+2);

                for (uint64_t i=0, offset=0, j=0; i<m_size; ++i, ++j) {
                    if (j==m_max_symbol) {
                        ++offset;
                        j = 0;
                    }
                    ++tmp[offset+input[i]*m_chunks];
                }

                for (uint64_t i=0,l=1; i<tmp.size(); ++i,++l) {
                    for (uint64_t j=0; j<tmp[i]; ++j) b[l++]=1;
                }

                //calc m_sigma
                bool write = true;
                uint64_t blocks=0;
                for (uint64_t i=1; i<b.size(); ++i) {
                    if (blocks==m_chunks) {
                        blocks = 0;
                        write = true;
                    }
                    if (b[i]) {
                        if (write) {
                            ++m_sigma;
                            write = false;
                        }
                    } else ++blocks;
                }

                m_bv_blocks = t_bitvector(std::move(b));
                util::init_support(m_bv_blocks_select1, &m_bv_blocks);
                util::init_support(m_bv_blocks_select0, &m_bv_blocks);
            }
            //calc inverse m_perm
            {
                bit_vector ipi(m_size,0);
                bit_vector pitmp(m_max_symbol,0);
                //calc m_bv_shortcut
                for (uint64_t i=0; i<m_chunks; ++i) {
                    for (uint64_t j = i*m_max_symbol,k=0; j<(i+1)*m_max_symbol and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_max_symbol,cycle_length=0;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    ipi[pos_pi]=1;
                                    steps = 0;
                                }
                                pos_pitmp = m_perm[pos_pi];
                                pos_pi = offset + m_perm[pos_pi];
                                ++steps;
                                ++cycle_length;
                            } while (!pitmp[pos_pitmp]);

                            if (steps==m_t and cycle_length>m_t) {
                                ipi[pos_pi]=1;
                            }
                        }
                    }
                    util::_set_zero_bits(pitmp);
                }
                m_bv_shortcut = t_bitvector(std::move(ipi));
                util::init_support(m_bv_shortcut_rank, &m_bv_shortcut);

                //calc m_shortcut
                m_shortcut = int_vector<>(m_bv_shortcut_rank(m_size),0,bits::hi(m_max_symbol-1)+1);

                for (uint64_t i=0; i<m_chunks; ++i) {
                    for (uint64_t j = i*m_max_symbol,k=0; j<(i+1)*m_max_symbol and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_max_symbol,cycle_length=0,back_pointer=k;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    m_shortcut[m_bv_shortcut_rank(pos_pi)]=back_pointer;
                                    back_pointer=pos_pitmp;
                                    steps = 0;
                                }
                                pos_pitmp = m_perm[pos_pi];
                                pos_pi = offset + m_perm[pos_pi];
                                ++steps;
                                ++cycle_length;
                            } while (!pitmp[pos_pitmp]);

                            if (steps==m_t and cycle_length>m_t) {
                                m_shortcut[m_bv_shortcut_rank(pos_pi)]=back_pointer;
                            }
                        }
                    }
                    util::_set_zero_bits(pitmp);
                }
            }
        }

        //! Swap operator
        void swap(wt_gmr_2& fs) {
            if (this != &fs) {
                m_bv_blocks.swap(fs.m_bv_blocks);
                m_bv_chunks.swap(fs.m_bv_chunks);
                m_bv_shortcut.swap(fs.m_bv_shortcut);
                m_perm.swap(fs.m_perm);
                m_shortcut.swap(fs.m_shortcut);
                util::swap_support(m_bv_blocks_select0, fs.m_bv_blocks_select0, &m_bv_blocks, &(fs.m_bv_blocks));
                util::swap_support(m_bv_blocks_select1, fs.m_bv_blocks_select1, &m_bv_blocks, &(fs.m_bv_blocks));
                util::swap_support(m_bv_chunks_select1, fs.m_bv_chunks_select1, &m_bv_chunks, &(fs.m_bv_chunks));
                util::swap_support(m_bv_chunks_select0, fs.m_bv_chunks_select0, &m_bv_chunks, &(fs.m_bv_chunks));
                util::swap_support(m_bv_shortcut_rank, fs.m_bv_shortcut_rank, &m_bv_shortcut, &(fs.m_bv_shortcut));
                std::swap(m_size, fs.m_size);
                std::swap(m_t, fs.m_t);
                std::swap(m_max_symbol,  fs.m_max_symbol);
                std::swap(m_chunks,  fs.m_chunks);
                std::swap(m_sigma,  fs.m_sigma);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        value_type operator[](size_type i)const {
            uint64_t chunk = i/m_max_symbol;
            bool jump=true;
            uint64_t x=i, value=i-(chunk*m_max_symbol);
            while (m_perm[x]!=value) {
                if (jump and m_bv_shortcut[x]==1) {
                    x  = m_shortcut[m_bv_shortcut_rank(x)]+(chunk*m_max_symbol);
                    jump = false;
                } else {
                    x = m_perm[x]+(chunk*m_max_symbol);
                }
            }
            return m_bv_chunks_select1(x+1)-x-(chunk*m_max_symbol)-1;
        }

        std::pair<size_type, value_type>	inverse_select(size_type i)const {
            uint64_t chunk = i/m_max_symbol;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_max_symbol);
            while (m_perm[x]!=value) {
                if (jump and m_bv_shortcut[x]==1) {
                    x  = m_shortcut[m_bv_shortcut_rank(x)]+(chunk*m_max_symbol);
                    jump = false;
                } else {
                    x = m_perm[x]+(chunk*m_max_symbol);
                }
            }
            uint64_t tmp = m_bv_chunks_select1(x+1);
            uint64_t c = tmp-x-(chunk*m_max_symbol)-1;

            uint64_t ones_before_c = m_bv_blocks_select0(c*m_chunks+1)-(c*m_chunks+1)+1;
            uint64_t c_before_chunk = m_bv_blocks_select0(c*m_chunks+chunk+1)-(c*m_chunks+chunk+1)+1-ones_before_c;
            uint64_t c_in_chunk = tmp-m_bv_chunks_select0(c+1+chunk*m_max_symbol)-1;
            return std::make_pair(c_before_chunk+c_in_chunk,c);
        }

        size_type select(size_type i, value_type c)const {

	    uint64_t ones_before_c = m_bv_blocks_select0(c*m_chunks+1)-(c*m_chunks);
            uint64_t chunk = m_bv_blocks_select1(ones_before_c+i)-ones_before_c-(c*m_chunks+1)-i+1;
            uint64_t c_ones_before_chunk = m_bv_blocks_select0(c*m_chunks+chunk+1)-(c*m_chunks+chunk)-ones_before_c;
            uint64_t pi_pos = m_bv_chunks_select0(chunk*m_max_symbol+c+1)+(i-c_ones_before_chunk)-chunk*m_max_symbol-c-1;

            return m_perm[pi_pos]+chunk*m_max_symbol;
        }

        size_type rank(size_type i, value_type c)const {

            if (c>m_max_symbol-1) return 0;
            if (i<=0) return 0;

            uint64_t chunk = (i-1)/m_max_symbol;
            uint64_t ones_before_c = m_bv_blocks_select0(c*m_chunks+1)-(c*m_chunks+1)+1;
            uint64_t c_ones_before_chunk = m_bv_blocks_select0(c*m_chunks+chunk+1)-(c*m_chunks+chunk+1)+1-ones_before_c;

            uint64_t c_ones_in_chunk = 0;

            size_type search_begin = m_bv_chunks_select0(chunk*m_max_symbol+1+c)-(chunk*m_max_symbol+1+c)+1;
            size_type search_end = m_bv_chunks_select0(chunk*m_max_symbol+2+c)-(chunk*m_max_symbol+2+c)+1;

            size_type val = (i-1)%m_max_symbol;
            if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                while (search_begin < search_end and m_perm[search_begin] <= val) {
                    ++search_begin;
                    ++c_ones_in_chunk;
                }
            } else {
                c_ones_in_chunk = lower_bound(m_perm.begin()+search_begin, m_perm.begin()+search_end, val+1)-m_perm.begin()-search_begin;
            }
            return c_ones_before_chunk+c_ones_in_chunk;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_max_symbol, out, child, "sigma");
            written_bytes += write_member(m_chunks, out, child, "chunks");
            written_bytes += write_member(m_t, out, child, "t");
            written_bytes += write_member(m_sigma, out, child, "test_sigma");
            written_bytes += m_bv_blocks.serialize(out, child, "m_bv_blocks");
            written_bytes += m_bv_chunks.serialize(out, child, "bv_chunks");
            written_bytes += m_bv_shortcut.serialize(out, child, "bv_shortcut");
            written_bytes += m_perm.serialize(out, child, "m_perm");
            written_bytes += m_shortcut.serialize(out, child, "m_shortcut");
            written_bytes += m_bv_blocks_select0.serialize(out, child, "m_bv_blocks_select0");
            written_bytes += m_bv_blocks_select1.serialize(out, child, "m_bv_blocks_select1");
            written_bytes += m_bv_chunks_select0.serialize(out, child, "m_bv_chunks_select0");
            written_bytes += m_bv_chunks_select1.serialize(out, child, "m_bv_chunks_select1");
            written_bytes += m_bv_shortcut_rank.serialize(out, child, "m_bv_shortcut_rank");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_max_symbol, in);
            read_member(m_chunks, in);
            read_member(m_t, in);
            read_member(m_sigma, in);
            m_bv_blocks.load(in);
            m_bv_chunks.load(in);
            m_bv_shortcut.load(in);
            m_perm.load(in);
            m_shortcut.load(in);
            m_bv_blocks_select0.load(in, &m_bv_blocks);
            m_bv_blocks_select1.load(in, &m_bv_blocks);
            m_bv_chunks_select0.load(in, &m_bv_chunks);
            m_bv_chunks_select1.load(in, &m_bv_chunks);
            m_bv_shortcut_rank.load(in, &m_bv_shortcut);
        }
};

template<class t_bitvector = bit_vector,
         class t_select = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_rank = typename t_bitvector::rank_1_type>
class wt_gmr_3
{
    private:

        t_bitvector m_bv_blocks; // 0 indicates end of block. Corresponds to B in the paper.
        t_bitvector m_bv_chunks; // 0 indicates end of symbol in chunk. Corresponds to X in the paper.
        t_bitvector m_bv_shortcut; // 1 indicates a shortcut

        int_vector<> m_shortcut; // Contains the shortcuts. Corresponds to S in the paper.
        int_vector<> m_perm; // Contains permutation of each chunk. Corresponds to \f$ \pi \f$ in the paper.

        t_rank m_bv_shortcut_rank;
        t_select m_bv_blocks_select1, m_bv_chunks_select1;
        t_select_zero m_bv_blocks_select0, m_bv_chunks_select0;

        uint64_t m_size; // input length
        uint64_t m_t = 32; // shortcut value
        uint64_t m_max_symbol = 0; // maximum character + 1
        uint64_t m_chunks; // number of chunks
        uint64_t m_chunksize;
        uint64_t m_sigma = 0 ;

    public:
        typedef int_vector<>::size_type size_type;
        typedef int_vector<>::value_type value_type;
        typedef wt_tag index_category;
        typedef int_alphabet_tag alphabet_category;
        enum {lex_ordered=0};

        const size_type&       sigma = m_sigma; // Todo

        wt_gmr_3() {}

        template<uint8_t int_width>
        wt_gmr_3(int_vector_buffer<int_width>& input, size_type size) : m_size(size) {
            // Determin max. symbol
            for (uint64_t i=0; i<m_size; ++i) {
                if (m_max_symbol < input[i]) m_max_symbol = input[i];
            }
            ++m_max_symbol;

            m_chunksize = (1 << (bits::hi(m_max_symbol-1)+1));
            m_chunks = (m_size+m_chunksize-1)/m_chunksize;

            m_perm = int_vector<>(m_size,0,bits::hi(m_max_symbol-1)+1);
            {
                uint64_t x_pos = 0;
                bit_vector x(m_size+m_chunks*m_max_symbol+1,0);

                //fill m_perm and m_bv_chunks for every chunk
                for (uint64_t i=0; i<m_chunks; ++i) {
                    int_vector<> symbols(m_max_symbol,0,bits::hi(m_max_symbol-1)+2);

                    //calc symbols
                    for (uint64_t j=i*m_chunksize; j<(i+1)*m_chunksize and j<m_size; ++j) {
                        ++symbols[input[j]];
                    }
                    //calc m_bv_chunks
                    for (uint64_t j=0; j<m_max_symbol; ++j,++x_pos) {
                        for (uint64_t k=0; k<symbols[j]; ++k) x[++x_pos]=1;
                    }

                    //calc symbols prefix sum
                    for (uint64_t j=0, tmp=0, sum=0; j<m_max_symbol; ++j) {
                        tmp = symbols[j];
                        symbols[j] = sum;
                        sum += tmp;
                    }
                    //calc m_perm
                    for (uint64_t j=i*m_chunksize, k=0; j<(i+1)*m_chunksize and j<m_size; ++j,++k) {
                        m_perm[i*m_chunksize+(symbols[input[j]]++)] = k;
                    }
                }
                m_bv_chunks = t_bitvector(std::move(x));
                util::init_support(m_bv_chunks_select1, &m_bv_chunks);
                util::init_support(m_bv_chunks_select0, &m_bv_chunks);

            }
            //calc m_bv_blocks
            {
                bit_vector b(m_size+m_max_symbol*m_chunks+1,0);
                int_vector<> tmp(m_max_symbol*m_chunks,0,bits::hi(m_max_symbol-1)+2);

                for (uint64_t i=0, offset=0, j=0; i<m_size; ++i, ++j) {
                    if (j==m_chunksize) {
                        ++offset;
                        j = 0;
                    }
                    ++tmp[offset+input[i]*m_chunks];
                }

                for (uint64_t i=0,l=1; i<tmp.size(); ++i,++l) {
                    for (uint64_t j=0; j<tmp[i]; ++j) b[l++]=1;
                }

                //calc m_sigma
                bool write = true;
                uint64_t blocks=0;
                for (uint64_t i=1; i<b.size(); ++i) {
                    if (blocks==m_chunks) {
                        blocks = 0;
                        write = true;
                    }
                    if (b[i]) {
                        if (write) {
                            ++m_sigma;
                            write = false;
                        }
                    } else ++blocks;
                }
                m_bv_blocks = t_bitvector(std::move(b));
                util::init_support(m_bv_blocks_select1, &m_bv_blocks);
                util::init_support(m_bv_blocks_select0, &m_bv_blocks);
            }

            //calc inverse m_perm
            {
                bit_vector ipi(m_size,0);
                bit_vector pitmp(m_chunksize,0);
                //calc m_bv_shortcut
                for (uint64_t i=0; i<m_chunks; ++i) {
                    for (uint64_t j = i*m_chunksize,k=0; j<(i+1)*m_chunksize and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_chunksize,cycle_length=0;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    ipi[pos_pi]=1;
                                    steps = 0;
                                }
                                pos_pitmp = m_perm[pos_pi];
                                pos_pi = offset + m_perm[pos_pi];
                                ++steps;
                                ++cycle_length;
                            } while (!pitmp[pos_pitmp]);

                            if (steps==m_t and cycle_length>m_t) {
                                ipi[pos_pi]=1;
                            }
                        }
                    }
                    util::_set_zero_bits(pitmp);
                }
                m_bv_shortcut = t_bitvector(std::move(ipi));
                util::init_support(m_bv_shortcut_rank, &m_bv_shortcut);

                //calc m_shortcut
                m_shortcut= int_vector<>(m_bv_shortcut_rank(m_size),0,bits::hi(m_max_symbol-1)+1);

                for (uint64_t i=0; i<m_chunks; ++i) {
                    for (uint64_t j = i*m_chunksize,k=0; j<(i+1)*m_chunksize and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_chunksize,cycle_length=0,back_pointer=k;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    m_shortcut[m_bv_shortcut_rank(pos_pi)]=back_pointer;
                                    back_pointer=pos_pitmp;
                                    steps = 0;
                                }
                                pos_pitmp = m_perm[pos_pi];
                                pos_pi = offset + m_perm[pos_pi];
                                ++steps;
                                ++cycle_length;
                            } while (!pitmp[pos_pitmp]);

                            if (steps==m_t and cycle_length>m_t) {
                                m_shortcut[m_bv_shortcut_rank(pos_pi)]=back_pointer;
                            }
                        }
                    }
                    util::_set_zero_bits(pitmp);
                }
            }
        }

        //! Swap operator
        void swap(wt_gmr_3& fs) {
            if (this != &fs) {
                m_bv_blocks.swap(fs.m_bv_blocks);
                m_bv_chunks.swap(fs.m_bv_chunks);
                m_bv_shortcut.swap(fs.m_bv_shortcut);
                m_perm.swap(fs.m_perm);
                m_shortcut.swap(fs.m_shortcut);
                util::swap_support(m_bv_blocks_select0, fs.m_bv_blocks_select0, &m_bv_blocks, &(fs.m_bv_blocks));
                util::swap_support(m_bv_blocks_select1, fs.m_bv_blocks_select1, &m_bv_blocks, &(fs.m_bv_blocks));
                util::swap_support(m_bv_chunks_select1, fs.m_bv_chunks_select1, &m_bv_chunks, &(fs.m_bv_chunks));
                util::swap_support(m_bv_chunks_select0, fs.m_bv_chunks_select0, &m_bv_chunks, &(fs.m_bv_chunks));
                util::swap_support(m_bv_shortcut_rank, fs.m_bv_shortcut_rank, &m_bv_shortcut, &(fs.m_bv_shortcut));
                std::swap(m_size, fs.m_size);
                std::swap(m_t, fs.m_t);
                std::swap(m_max_symbol,  fs.m_max_symbol);
                std::swap(m_chunks,  fs.m_chunks);
                std::swap(m_chunksize,  fs.m_chunksize);
                std::swap(m_sigma,  fs.m_sigma);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        value_type operator[](size_type i)const {
            uint64_t chunk = i/m_chunksize;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_chunksize);
            while (m_perm[x]!=value) {
                if (jump and m_bv_shortcut[x]==1) {
                    x  = m_shortcut[m_bv_shortcut_rank(x)]+(chunk*m_chunksize);
                    jump = false;
                } else {
                    x = m_perm[x]+(chunk*m_chunksize);
                }
            }
            return m_bv_chunks_select1(x+1)-x-(chunk*m_max_symbol)-1;
        }

        std::pair<size_type, value_type>	inverse_select(size_type i)const {
            uint64_t chunk = i/m_chunksize;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_chunksize);
            while (m_perm[x]!=value) {
                if (jump and m_bv_shortcut[x]==1) {
                    x  = m_shortcut[m_bv_shortcut_rank(x)]+(chunk*m_chunksize);
                    jump = false;
                } else {
                    x = m_perm[x]+(chunk*m_chunksize);
                }
            }
            uint64_t tmp = m_bv_chunks_select1(x+1);
            uint64_t c = tmp-x-(chunk*m_max_symbol)-1;

	    uint64_t ones_before_c = m_bv_blocks_select0(c*m_chunks+1)-(c*m_chunks+1)+1;
	    uint64_t c_before_chunk = m_bv_blocks_select0(c*m_chunks+chunk+1)-(c*m_chunks+chunk+1)+1-ones_before_c;
            uint64_t c_in_chunk = tmp-m_bv_chunks_select0(c+1+chunk*m_max_symbol)-1;
            return std::make_pair(c_before_chunk+c_in_chunk,c);
        }

        size_type select(size_type i, value_type c)const {

	    uint64_t ones_before_c = m_bv_blocks_select0(c*m_chunks+1)-(c*m_chunks);
	    uint64_t chunk = m_bv_blocks_select1(ones_before_c+i)-ones_before_c-(c*m_chunks+1)-i+1;
	    uint64_t c_ones_before_chunk = m_bv_blocks_select0(c*m_chunks+chunk+1)-(c*m_chunks+chunk)-ones_before_c;
	    uint64_t pi_pos = m_bv_chunks_select0(chunk*m_max_symbol+c+1)+(i-c_ones_before_chunk)-chunk*m_max_symbol-c-1;

            return m_perm[pi_pos]+chunk*m_chunksize;
        }

        size_type rank(size_type i, value_type c)const {

            if (c>m_max_symbol-1) return 0;
            if (i<=0) return 0;

            uint64_t chunk = (i-1)/m_chunksize;
            uint64_t ones_before_c = m_bv_blocks_select0(c*m_chunks+1)-(c*m_chunks+1)+1;
            uint64_t c_ones_before_chunk = m_bv_blocks_select0(c*m_chunks+chunk+1)-(c*m_chunks+chunk+1)+1-ones_before_c;

            uint64_t c_ones_in_chunk = 0;

            size_type search_begin = m_bv_chunks_select0(chunk*m_max_symbol+1+c)-(chunk*m_max_symbol+1+c)+1;
            size_type search_end = m_bv_chunks_select0(chunk*m_max_symbol+2+c)-(chunk*m_max_symbol+2+c)+1;

            size_type val = (i-1)%m_chunksize;
            if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                while (search_begin < search_end and m_perm[search_begin] <= val) {
                    ++search_begin;
                    ++c_ones_in_chunk;
                }
            } else {
                c_ones_in_chunk = lower_bound(m_perm.begin()+search_begin, m_perm.begin()+search_end, val+1)-m_perm.begin()-search_begin;
            }
            return c_ones_before_chunk+c_ones_in_chunk;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_max_symbol, out, child, "sigma");
            written_bytes += write_member(m_chunks, out, child, "chunks");
            written_bytes += write_member(m_chunksize, out, child, "chunksize");
            written_bytes += write_member(m_t, out, child, "t");
            written_bytes += write_member(m_sigma, out, child, "test_sigma");
            written_bytes += m_bv_blocks.serialize(out, child, "m_bv_blocks");
            written_bytes += m_bv_chunks.serialize(out, child, "bv_chunks");
            written_bytes += m_bv_shortcut.serialize(out, child, "bv_shortcut");
            written_bytes += m_perm.serialize(out, child, "m_perm");
            written_bytes += m_shortcut.serialize(out, child, "m_shortcut");
            written_bytes += m_bv_blocks_select0.serialize(out, child, "m_bv_blocks_select0");
            written_bytes += m_bv_blocks_select1.serialize(out, child, "m_bv_blocks_select1");
            written_bytes += m_bv_chunks_select0.serialize(out, child, "m_bv_chunks_select0");
            written_bytes += m_bv_chunks_select1.serialize(out, child, "m_bv_chunks_select1");
            written_bytes += m_bv_shortcut_rank.serialize(out, child, "m_bv_shortcut_rank");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_max_symbol, in);
            read_member(m_chunks, in);
            read_member(m_chunksize, in);
            read_member(m_t, in);
            read_member(m_sigma, in);
            m_bv_blocks.load(in);
            m_bv_chunks.load(in);
            m_bv_shortcut.load(in);
            m_perm.load(in);
            m_shortcut.load(in);
            m_bv_blocks_select0.load(in, &m_bv_blocks);
            m_bv_blocks_select1.load(in, &m_bv_blocks);
            m_bv_chunks_select0.load(in, &m_bv_chunks);
            m_bv_chunks_select1.load(in, &m_bv_chunks);
            m_bv_shortcut_rank.load(in, &m_bv_shortcut);
        }
};

}

#endif
