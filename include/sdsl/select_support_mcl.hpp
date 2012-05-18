/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file select_support_mcl.hpp
    \brief select_support_mcl.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_JMC
#define INCLUDED_SDSL_SELECT_SUPPORT_JMC

#include "int_vector.hpp"
#include "select_support.hpp"

//#define SDSL_DEBUG_SELECT_SUPPORT_JMC

#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
#include "testutils.hpp"
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<uint8_t bit_pattern, uint8_t pattern_len>
struct select_support_mcl_trait {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return 0;
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return 0;
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return 0;
    }

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return 0;
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return 0;
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        return 0;
    }

    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return 0;
    }
    static uint64_t get_carry(uint64_t w) {
        return 0;
    }
};

template<>
struct select_support_mcl_trait<0,1> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v) {
        return v.bit_size()-util::get_one_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bit_magic::b1Cnt((~w)& bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bit_magic::i1BP(~w & bit_magic::Li0Mask[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bit_magic::b1Cnt(~w);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bit_magic::i1BP(~w, i);
    }
    static bool found_arg(size_type i, const bit_vector& v) {
        return !v[i];
    }
    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return 0;
    }
    static uint64_t get_carry(uint64_t w) {
        return 0;
    }
};

template<>
struct select_support_mcl_trait<1,1> {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return util::get_one_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bit_magic::b1Cnt(w & bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bit_magic::i1BP(w & bit_magic::Li0Mask[offset], i);
    }

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bit_magic::b1Cnt(w);
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bit_magic::i1BP(w, i);
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        return v[i];
    }

    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return 0;
    }
    static uint64_t get_carry(uint64_t w) {
        return 0;
    }
};

template<>
struct select_support_mcl_trait<10,2> {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return util::get_onezero_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bit_magic::b1Cnt(bit_magic::b10Map(w, carry) & bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b10Map(w, carry) & bit_magic::Li0Mask[offset], i);
    }

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bit_magic::b10Cnt(w, carry);
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b10Map(w, carry), i);
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        if (i > 0 and v[i-1] and !v[i])
            return true;
        return false;
    }

    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return word_pos ? (*(data-1)>>63) : 0;
    }
    static uint64_t get_carry(uint64_t w) {
        return w>>63;
    }
};

template<>
struct select_support_mcl_trait<01,2> {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return util::get_zeroone_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bit_magic::b1Cnt(bit_magic::b01Map(w, carry) & bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b01Map(w, carry) & bit_magic::Li0Mask[offset], i);
    }

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bit_magic::b01Cnt(w, carry);
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b01Map(w, carry), i);
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        if (i > 0 and !v[i-1] and v[i])
            return true;
        return false;
    }

    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return word_pos ? (*(data-1)>>63) : 1;
    }
    static uint64_t get_carry(uint64_t w) {
        return w>>63;
    }
};


//! A class supporting constant time select queries (proposed by Munro/Clark, 1996) enhanced by broadword computing tricks.
/*!
 * \par Space usage
 *      The space usage of the data structure depends on the number of \f$ m \f$ of ones in the
 *      original bitvector $b$. We store the position of every $4096$th set bit
 *      (called L1-sampled bits) of $b$.
 *      This takes in the worst case \f$\frac{m}{4096} \log{n} \leq \frac{64}{n}\f$ bits.
 *      Next,
 *      (1) if the distance of two adjacent L1-sampled bits $b[i]$ and $b[j]$
 *      is greater or equal than $\log^4 n$, then
 *      we store each of the 4096 positions of the set $b$ in [i..j-1] with
 *      $\log{n}$ bits. This results in at most
 *      \$ \frac{4096\cdot \log n}{\log^4 n}=\frac{4096}{\log^3 n}\$ bits per bit.
 *      For a bitvector of 1GB, i.e. \f$ \log n = 35 \f$ we get about 0.01 bits per bit.
 *      If the $j-i+1 < \log^4 n$ then
 *      (2) we store the relative position of every $64$th set bit (called L2-sampled bits)
 *      in b[i..j-1] in at most $4\log\log n$ bits per L2-sampled bits.
 *      An pessimistic upper bound for the space would be
 *      \f$ \frac{4\log\log n}{64} \leq \frac{24}{64} = 0.375\f$ bit per
 *      bit (since $\log\log n\leq 6$. It is very pessimistic, since we store
 *      the relative position in $\log\log(j-i+1)\leq \log\log n$ bits.
 * @ingroup select_support_group
 */
template<uint8_t b=1, uint8_t pattern_len=1>
class select_support_mcl : public select_support
{
    public:
        typedef bit_vector bit_vector_type;
    private:
        uint32_t m_logn, m_logn2, m_logn4; // log(size of the supported bit_vector), logn2=\f$ (\log n)^2\f$, \f$logn4=(\log n)^4\f$
        int_vector<0> m_superblock;        // there exists \f$\frac{n}{4096} < \frac{n}{\log n}\f$ entries
        // entry i equals the answer to select_1(B,i*4096)
        int_vector<0>* m_longsuperblock;
        int_vector<0>* m_miniblock;
        size_type m_arg_cnt;
        void copy(const select_support_mcl<b, pattern_len>& ss);
        void construct();
        void initData();
        void init_fast(const int_vector<1>* v=NULL);
    public:
        select_support_mcl(const int_vector<1>* v=NULL);
        select_support_mcl(const select_support_mcl<b,pattern_len>& ss);
        ~select_support_mcl();
        void init(const int_vector<1>* v=NULL);
        void init_slow(const int_vector<1>* v=NULL);
        //! Select function
        /*! \sa select_support.select
         */
        inline const size_type select(size_type i) const;
        //! Alias for select(i).
        inline const size_type operator()(size_type i)const;
        size_type serialize(std::ostream& out)const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v=NULL);
        select_support_mcl<b, pattern_len>& operator=(const select_support_mcl& ss);

        //! Swap operator
        /*! This swap Operator swaps two select_support_mcls in constant time.
         *  All members (excluded the pointer to the supported SDSBitVector) are swapped.
         */
        void swap(select_support_mcl<b, pattern_len>& ss);
        //! Equality Operator
        /*! Two select_support_mcls are equal if all member variables are equal.
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator!=
         */
        bool operator==(const select_support_mcl<b, pattern_len>& ss)const;
        //! Unequality Operator
        /*! Two select_support_mcls are not equal if any member variable are not equal.
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator==
         */
        bool operator!=(const select_support_mcl<b, pattern_len>& ss)const;
#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label="select";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<")\n";
        }
#endif
};


template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::construct()
{
    m_longsuperblock	= NULL;
    m_miniblock			= NULL;
    m_arg_cnt			= 0;
}

template<uint8_t b, uint8_t pattern_len>
select_support_mcl<b,pattern_len>::select_support_mcl(const int_vector<1>* f_v):select_support(f_v)
{
    construct();
    init(f_v);
}

template<uint8_t b, uint8_t pattern_len>
select_support_mcl<b,pattern_len>::select_support_mcl(const select_support_mcl& ss):select_support(ss.m_v)
{
    construct();
    copy(ss);
}

template<uint8_t b, uint8_t pattern_len>
select_support_mcl<b, pattern_len>& select_support_mcl<b,pattern_len>::operator=(const select_support_mcl& ss)
{
    if (this != &ss) {
        copy(ss);
    }
    return *this;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::swap(select_support_mcl& ss)
{
    std::swap(m_logn, ss.m_logn);
    std::swap(m_logn2, ss.m_logn2);
    std::swap(m_logn4, ss.m_logn4);
    m_superblock.swap(ss.m_superblock);
    std::swap(m_longsuperblock, ss.m_longsuperblock);
    std::swap(m_miniblock, ss.m_miniblock);
    std::swap(m_arg_cnt, ss.m_arg_cnt);
//	std::swap(m_v, ss.m_v); // see rank_support_v swap.
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::copy(const select_support_mcl<b, pattern_len>& ss)
{
    m_logn			= ss.m_logn;		// copy log n
    m_logn2			= ss.m_logn2;		// copy (logn)^2
    m_logn4			= ss.m_logn4;		// copy (logn)^4
    m_superblock	= ss.m_superblock;  // copy long superblock
    m_arg_cnt		= ss.m_arg_cnt;	    // copy count of 1-bits
    m_v				= ss.m_v;		    // copy pointer to the supported bit vector
    size_type sb	= (m_arg_cnt+4095)>>12;
    if (m_longsuperblock!=NULL)
        delete [] m_longsuperblock;
    m_longsuperblock = NULL;
    if (ss.m_longsuperblock!=NULL) {
        m_longsuperblock = new int_vector<0>[sb]; //copy longsuperblocks
        for (size_type i=0; i<sb; ++i) {
            m_longsuperblock[i] = ss.m_longsuperblock[i];
        }
    }
    if (m_miniblock!=NULL)
        delete [] m_miniblock;
    m_miniblock = NULL;
    if (ss.m_miniblock!=NULL) {
        m_miniblock = new int_vector<0>[sb]; // copy miniblocks
        for (size_type i=0; i<sb; ++i) {
            m_miniblock[i] = ss.m_miniblock[i];
        }
    }
}

template<uint8_t b, uint8_t pattern_len>
select_support_mcl<b,pattern_len>::~select_support_mcl()
{
    if (m_longsuperblock!=NULL)
        delete[] m_longsuperblock;
    if (m_miniblock!=NULL)
        delete[] m_miniblock;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::init_slow(const int_vector<1>* v)
{
    set_vector(v);
    initData();
    if (m_v==NULL)
        return;
#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    stop_watch sw;
    sw.start();
#endif
    // Count the number of arguments in the bit vector
    m_arg_cnt = select_support_mcl_trait<b,pattern_len>::arg_cnt(*v);
#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    sw.stop();
    std::cerr<<"Count arguments takes "<<sw.get_real_time()<<" ms"<<std::endl;
    sw.start();
#endif

    const size_type SUPER_BLOCK_SIZE = 4096;

    if (m_arg_cnt==0) // if there are no arguments in the vector we are done...
        return;

    size_type sb = (m_arg_cnt+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; // number of superblocks
    if (m_miniblock != NULL) delete [] m_miniblock;
    m_miniblock = new int_vector<0>[sb];

    m_superblock = int_vector<0>(sb, 0, m_logn);// TODO: hier koennte man logn noch optimieren...s


    size_type arg_position[SUPER_BLOCK_SIZE], arg_cnt=0;
    size_type sb_cnt=0;
    for (size_type i=0; i < v->size(); ++i) {
        if (select_support_mcl_trait<b,pattern_len>::found_arg(i, *v)) {
            arg_position[ arg_cnt%SUPER_BLOCK_SIZE ] = i;
            assert(arg_position[arg_cnt%SUPER_BLOCK_SIZE] == i);
            ++arg_cnt;
            if (arg_cnt % SUPER_BLOCK_SIZE == 0 or arg_cnt == m_arg_cnt) { //
                assert(sb_cnt < sb);
                m_superblock[sb_cnt] = arg_position[0];

                size_type pos_diff = arg_position[(arg_cnt-1)%SUPER_BLOCK_SIZE]-arg_position[0];
                if (pos_diff > m_logn4) { // longblock
                    if (m_longsuperblock == NULL) m_longsuperblock = new int_vector<0>[sb]; // create longsuperblock
                    m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bit_magic::l1BP(arg_position[(arg_cnt-1)%SUPER_BLOCK_SIZE]) + 1);

                    for (size_type j=0; j <= (arg_cnt-1)%SUPER_BLOCK_SIZE ; ++j) m_longsuperblock[sb_cnt][j] = arg_position[j]; // copy argument positions to longsuperblock
                } else { // short block
                    m_miniblock[sb_cnt] = int_vector<0>(64, 0, bit_magic::l1BP(pos_diff)+1);
                    for (size_type j=0; j <= (arg_cnt-1)%SUPER_BLOCK_SIZE; j+=64) {
                        m_miniblock[sb_cnt][j/64] = arg_position[j]-arg_position[0];
                    }
                }
                ++sb_cnt;
            }
        }
    }
//	std::cout<<"m_superblock[0]="<<m_superblock[0]<<std::endl;
//	std::cout<<"sb_cnt="<<sb_cnt<<std::endl;
//	std::cout<<"arg_position[0]="<<arg_position[0]<<std::endl;
#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    sw.stop();
    std::cerr<<"Init tables  takes "<<sw.get_real_time()<<" ms"<<std::endl;
#endif
}

// TODO: find bug, detected by valgrind
template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::init_fast(const int_vector<1>* v)
{
    set_vector(v);
    initData();
    if (m_v==NULL)
        return;
//#define SDSL_DEBUG_SELECT_SUPPORT_JMC
#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    stop_watch sw;
    sw.start();
#endif
    // Count the number of arguments in the bit vector
    m_arg_cnt = select_support_mcl_trait<b,pattern_len>::arg_cnt(*v);
#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    sw.stop();
    std::cerr<<"Count arguments takes "<<sw.get_real_time()<<" ms"<<std::endl;
    std::cerr<<"m_arg_cnt="<<m_arg_cnt<<std::endl;
    sw.start();
#endif

    const size_type SUPER_BLOCK_SIZE = 64*64;

    if (m_arg_cnt==0) // if there are no arguments in the vector we are done...
        return;

//	size_type sb = (m_arg_cnt+63+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; // number of superblocks, add 63 as the last block could contain 63 uninitialized bits
    size_type sb = (m_arg_cnt+SUPER_BLOCK_SIZE-1)/SUPER_BLOCK_SIZE; // number of superblocks
    if (m_miniblock != NULL) delete [] m_miniblock;
    m_miniblock = new int_vector<0>[sb];

    m_superblock = int_vector<0>(sb, 0, m_logn);// TODO: hier koennte man logn noch optimieren...s

#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    std::cerr<<"H1"<<std::endl;
#endif

    bit_vector::size_type arg_position[SUPER_BLOCK_SIZE];
    const uint64_t* data = v->data();
    uint64_t carry_new=0;
    size_type last_k64 = 1, sb_cnt=0;
    for (size_type i=0, cnt_old=0, cnt_new=0, last_k64_sum=1; i < v->capacity(); i+=64, ++data) {
        cnt_new += select_support_mcl_trait<b, pattern_len>::args_in_the_word(*data, carry_new);
        if (cnt_new >= last_k64_sum) {
            arg_position[last_k64-1] = i + select_support_mcl_trait<b, pattern_len>::ith_arg_pos_in_the_word(*data, last_k64_sum  - cnt_old, carry_new);
            last_k64 += 64;
            last_k64_sum += 64;

            if (last_k64 == SUPER_BLOCK_SIZE+1) {
                m_superblock[sb_cnt] = arg_position[0];
                size_type pos_of_last_arg_in_the_block = arg_position[last_k64-65];

                for (size_type i=arg_position[last_k64-65]+1, j=last_k64-65; i < v->size() and j < SUPER_BLOCK_SIZE; ++i)
                    if (select_support_mcl_trait<b,pattern_len>::found_arg(i, *v)) {
                        pos_of_last_arg_in_the_block = i;
                        ++j;
                    }
                size_type pos_diff = pos_of_last_arg_in_the_block - arg_position[0];
                if (pos_diff > m_logn4) { // long block
                    if (m_longsuperblock == NULL) m_longsuperblock = new int_vector<0>[sb+1]; // create longsuperblock
                    // GEANDERT am 2010-07-17 +1 nach pos_of_last_arg..
                    m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bit_magic::l1BP(pos_of_last_arg_in_the_block) + 1);
                    for (size_type j=arg_position[0], k=0; j < pos_of_last_arg_in_the_block; ++j)
                        if (select_support_mcl_trait<b, pattern_len>::found_arg(j, *v)) {
//							if(k>=SUPER_BLOCK_SIZE){
//								std::cout<<"k="<<k<<" SUPER_BLOCK_SIZE="<<SUPER_BLOCK_SIZE<<std::endl;
//							}
                            m_longsuperblock[sb_cnt][k++] = j;
                        }
                } else {
                    m_miniblock[sb_cnt] = int_vector<0>(64, 0, bit_magic::l1BP(pos_diff)+1);
                    for (size_type j=0; j < SUPER_BLOCK_SIZE; j+=64) {
                        m_miniblock[sb_cnt][j/64] = arg_position[j]-arg_position[0];
                    }
                }
                ++sb_cnt;
                last_k64 = 1;
            }
        }
        cnt_old = cnt_new;
    }
    // TODO: handle last block!!! einfach ein longsuperblock anhaengen?
    if (last_k64 > 1) {
        // GEANDERT am 2010-09-27 +1 nach sb, da wir den superblock extra zu den existierenden anhaengen koennen
        if (m_longsuperblock == NULL) m_longsuperblock = new int_vector<0>[sb+1]; // create longsuperblock
        m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bit_magic::l1BP(v->size()-1) + 1);
        for (size_type i=arg_position[0],k=0; i < v->size(); ++i) {
            if (select_support_mcl_trait<b, pattern_len>::found_arg(i, *v)) {
                m_longsuperblock[sb_cnt][k++] = i;
            }
        }
        ++sb_cnt;
    }
#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
    sw.stop();
    std::cerr<<"Init tables  takes "<<sw.get_real_time()<<" ms"<<std::endl;
#endif
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::init(const int_vector<1>* v)
{
    if (pattern_len>1 or (v!=NULL and  v->size() < 100000))
        init_slow(v);
    else
        init_fast(v);
    return;
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_mcl<b,pattern_len>::size_type select_support_mcl<b,pattern_len>::select(size_type i)const
{
    assert(i > 0 and i <= m_arg_cnt);

    i = i-1;
    size_type sb_idx = i>>12;   // i/4096
    size_type offset = i&0xFFF; // i%4096
    if (m_longsuperblock!=NULL and !m_longsuperblock[sb_idx].empty()) {
//		std::cout<<"i="<<i<<" sb_idx="<<sb_idx<<" offset="<<offset<<std::endl;
//		std::cout<<"sb_size"<<m_longsuperblock[sb_idx].size()<<std::endl;
        return m_longsuperblock[sb_idx][offset];
    } else {
        if ((offset&0x3F)==0) {
            assert(sb_idx < m_superblock.size());
//			if( (offset>>6) >= m_miniblock[sb_idx].size() ){
//				std::cerr<<" i="<<i<<std::endl;
//				std::cerr<<" "<< (offset>>6) <<" >= "<<  m_miniblock[sb_idx].size() << std::endl;
//			}
            assert((offset>>6) < m_miniblock[sb_idx].size());
            return m_superblock[sb_idx] + m_miniblock[sb_idx][offset>>6/*/64*/];
        } else {
            i = i-(sb_idx<<12)-((offset>>6)<<6);
            // now i > 0 and i <= 64
            assert(i > 0);
//std::cout<<"sb_idx="<<sb_idx<<" "<<m_superblock.size()<<std::endl;
//std::cout<<"offset>>6="<< (offset>>6) << " "<<m_miniblock[sb_idx].size()<< std::endl;
            size_type pos = m_superblock[sb_idx] + m_miniblock[sb_idx][offset>>6] + 1;

            // now pos is the position from where we search for the ith argument
            size_type word_pos = pos>>6;
            size_type word_off = pos&0x3F;
//if( m_v->size() < word_pos*64 ){
//	std::cout<<"m_v->size()="<<m_v->size()<<" 64*word_pos="<<64*word_pos<<std::endl;
//}
            const uint64_t* data = m_v->data() + word_pos;
            uint64_t carry = select_support_mcl_trait<b,pattern_len>::init_carry(data, word_pos);
            size_type args = select_support_mcl_trait<b,pattern_len>::args_in_the_first_word(*data, word_off, carry);

            if (args >= i) {
                return (word_pos<<6)+select_support_mcl_trait<b,pattern_len>::ith_arg_pos_in_the_first_word(*data, i, word_off, carry);
            }
            word_pos+=1;
            size_type sum_args = args;
            carry = select_support_mcl_trait<b,pattern_len>::get_carry(*data);
            uint64_t old_carry = carry;
            args = select_support_mcl_trait<b,pattern_len>::args_in_the_word(*(++data), carry);
            while (sum_args + args < i) {
                sum_args += args;
                assert(data+1 < m_v->data() + (m_v->capacity()>>6));
                old_carry = carry;
                args = select_support_mcl_trait<b,pattern_len>::args_in_the_word(*(++data), carry);
                word_pos+=1;
            }
            return (word_pos<<6) +
                   select_support_mcl_trait<b,pattern_len>::ith_arg_pos_in_the_word(*data, i-sum_args, old_carry);
        }
    }
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_mcl<b,pattern_len>::size_type select_support_mcl<b,pattern_len>::operator()(size_type i)const
{
    return select(i);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::initData()
{
    m_arg_cnt = 0;
    if (m_v==NULL) {
        m_logn = m_logn2 = m_logn4 = 0;
    } else {
        m_logn = bit_magic::l1BP(m_v->capacity())+1; // TODO maybe it's better here to take a max(...,12)
        m_logn2 = m_logn*m_logn;
        m_logn4 = m_logn2*m_logn2;
    }
    if (m_longsuperblock!=NULL)
        delete[] m_longsuperblock;
    m_longsuperblock  = NULL;
    if (m_miniblock!=NULL)
        delete[] m_miniblock;
    m_miniblock 		= NULL;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::set_vector(const int_vector<1>* v)
{
    m_v = v;
}

template<uint8_t b, uint8_t pattern_len>
typename select_support_mcl<b,pattern_len>::size_type select_support_mcl<b,pattern_len>::serialize(std::ostream& out)const
{
    size_type written_bytes = 0;
    // write the number of 1-bits in the supported SDSBitVector
    out.write((char*) &m_arg_cnt, sizeof(size_type)/sizeof(char));
    written_bytes = sizeof(size_type)/sizeof(char);
    // number of superblocks in the data structure
    size_type sb = (m_arg_cnt+4095)>>12;

    if (m_arg_cnt) { // if there exists 1-bits to be supported
        written_bytes += m_superblock.serialize(out); // serialize superblocks
        int_vector<1> mini_or_long;// Helper vector: mini or long block?
        if (m_longsuperblock!=NULL) {
            mini_or_long.resize(sb); // resize indicatior BitVector to the number of superblocks
            for (size_type i=0; i< sb; ++i)
                mini_or_long[i] = !m_miniblock[i].empty();
        }
        written_bytes += mini_or_long.serialize(out);
        for (size_type i=0; i < sb; ++i)
            if (!mini_or_long.empty() and !mini_or_long[i]) {
                written_bytes += m_longsuperblock[i].serialize(out);
//		   			SDSCoder::compress(m_longsuperblock[i], z);
//					z.serialize(out);
//SDSBitVector *z = SDSCoder::compress(m_longsuperblock[i]);
//if(z!=NULL) delete z;
            } else {
                written_bytes += m_miniblock[i].serialize(out);
//		   			SDSCoder::compress(m_miniblock[i], z);
//		   			SDSCoder::decompress(z, zz);
//					assert(m_miniblock[i]==zz);
//					z.serialize(out);
//SDSBitVector *z = SDSCoder::compress(m_miniblock[i]);
//if(z!=NULL) delete z;
            }
    }
    return written_bytes;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_mcl<b,pattern_len>::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);
    initData();
    // read the number of 1-bits in the supported SDSBitVector
    in.read((char*) &m_arg_cnt, sizeof(size_type)/sizeof(char));
    size_type sb = (m_arg_cnt+4095)>>12;

    if (m_arg_cnt) { // if there exists 1-bits to be supported
//			SDSBitVector z;
//			z.load(in);// load compressed superblocks
//			SDSCoder::decompress(z, m_superblock);
        m_superblock.load(in); // load superblocks

        if (m_miniblock!=NULL) {
            delete[] m_miniblock;
            m_miniblock = NULL;
        }
        if (m_longsuperblock!=NULL) {
            delete[] m_longsuperblock;
            m_longsuperblock = NULL;
        }

        int_vector<1> mini_or_long;// Helper vector: mini or long block?
        mini_or_long.load(in); // Load the helper vector
        m_miniblock = new int_vector<0>[sb]; // Create miniblock int_vector<0>
        if (!mini_or_long.empty())
            m_longsuperblock = new int_vector<0>[sb]; // Create longsuperblock int_vector<0>

        for (size_type i=0; i < sb; ++i)
            if (!mini_or_long.empty() and not mini_or_long[i]) {
                m_longsuperblock[i].load(in);
//					z.load(in);
//					SDSCoder::decompress(z,m_longsuperblock[i]);
            } else {
                m_miniblock[i].load(in);
//					z.load(in);
//					SDSCoder::decompress(z,m_miniblock[i]);
            }
    }
}

template<uint8_t b, uint8_t pattern_len>
bool select_support_mcl<b,pattern_len>::operator==(const select_support_mcl& ss)const
{
    if (this == &ss)
        return true;
    if (m_logn == ss.m_logn and m_logn2 == ss.m_logn2 and m_logn4 == ss.m_logn4
        and m_arg_cnt == ss.m_arg_cnt and  m_superblock == ss.m_superblock) {
        if (m_longsuperblock == NULL and ss.m_longsuperblock != NULL)
            return false;
        if (ss.m_longsuperblock == NULL and m_longsuperblock != NULL)
            return false;
        if (m_miniblock == NULL and ss.m_miniblock != NULL)
            return false;
        if (ss.m_miniblock == NULL and m_miniblock != NULL)
            return false;
        size_type sb		= (m_arg_cnt+4095)>>12;
        if (m_miniblock != NULL)
            for (size_type i=0; i < sb; ++i)
                if (m_miniblock[i] != ss.m_miniblock[i])
                    return false;
        if (m_longsuperblock != NULL)
            for (size_type i=0; i < sb; ++i)
                if (m_longsuperblock[i] != ss.m_longsuperblock[i])
                    return false;
        if (m_v == NULL and ss.m_v != NULL)
            return false;
        if (m_v != NULL and ss.m_v == NULL)
            return false;
        if (m_v == NULL and ss.m_v == NULL)
            return false;
        return *ss.m_v == *m_v;// supported vectors have to be equal. See docu!!!
    } else
        return false;
}

template<uint8_t b, uint8_t pattern_len>
bool select_support_mcl<b,pattern_len>::operator!=(const select_support_mcl& ss)const
{
    return !(*this == ss);
}


}

#endif
