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
/*! \file csa_sada.hpp
    \brief csa_sada.hpp contains an implementation of the compressed suffix array.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_SADA
#define INCLUDED_SDSL_CSA_SADA

#include "enc_vector.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "suffixarrays.hpp"
#include "suffixarray_helper.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "bwt_construct.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens,  uint8_t fixedIntWidth>
class csa_sada;

template<uint8_t fixedIntWidth>
struct csa_sada_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct csa_sada_trait<32> {
    typedef int_vector<32> int_vector_type;
};

template<>
struct csa_sada_trait<64> {
    typedef int_vector<64> int_vector_type;
};


//! A class for the Compressed Suffix Array (CSA) proposed by Sadakane for practical implementation.
/*! The CSA is parameterized with an EncVector and the sample density SampleDens (\f$s_{SA}\f$).
  * I.e. every \f$s_{SA}th\f$ value from the original suffix array is explicitely stored with \f$\log n\f$ bits.
  *
  * The EncVector (default is sdsl::enc_vector) holds the \f$\Psi\f$-function and can be parametrized with \f$s_{\Psi}\f$.
  * \todo example, code example
  *  \sa csa_sada_theo
  * @ingroup csa
 */
template<class EncVector = enc_vector<>, uint32_t SampleDens = 32, uint32_t InvSampleDens = 64,  uint8_t fixedIntWidth = 0>
class csa_sada
{
    public:
        typedef uint64_t											      value_type;	// STL Container requirement
        typedef random_access_const_iterator<csa_sada> 				      const_iterator;// STL Container requirement
        typedef const_iterator 										      iterator;		// STL Container requirement
        typedef const value_type									      const_reference;
        typedef const_reference										      reference;
        typedef const_reference*									      pointer;
        typedef const pointer										      const_pointer;
        typedef int_vector<>::size_type								      size_type;		// STL Container requirement
        typedef size_type 											      csa_size_type;
        typedef ptrdiff_t  											      difference_type; // STL Container requirement
        typedef EncVector											      enc_vector_type;
        typedef psi_of_csa_psi<csa_sada>						 	      psi_type;
        typedef bwt_of_csa_psi<csa_sada>						 	      bwt_type;
        typedef const unsigned char*						 		      pattern_type;
        typedef unsigned char										      char_type;
        typedef typename csa_sada_trait<fixedIntWidth>::int_vector_type   sa_sample_type;
        typedef typename csa_sada_trait<fixedIntWidth>::int_vector_type   isa_sample_type;

        typedef csa_tag													  index_category;

        enum { sa_sample_dens = SampleDens,
               isa_sample_dens = InvSampleDens
             };

        friend class psi_of_csa_psi<csa_sada>;
        friend class bwt_of_csa_psi<csa_sada>;

//	static const uint32_t sample_dens = SampleDens;
    private:
        enc_vector_type m_psi;  // psi function
        psi_type m_psi_wrapper;
        bwt_type m_bwt;
        sa_sample_type 	m_sa_sample; // suffix array samples
        isa_sample_type m_isa_sample; // inverse suffix array samples
        int_vector<8>	m_char2comp; // =0 for the 0-byte and all characters which do not occur in the text
        int_vector<8> 	m_comp2char;
        int_vector<64>  m_C;
        uint16_t			m_sigma;

        uint64_t* m_psi_buf; // buffer for decoded psi values : TODO statisch uint64_t[groesse] anlegen, siehe: wavelet_tree fuer int_vector_file_buffer

        template<typename RandomAccessContainer>
        void construct(const RandomAccessContainer& sa, const unsigned char* str);

        template<typename RandomAccessContainer>
        void construct(RandomAccessContainer& sa, const unsigned char* str);

        template<typename RandomAccessContainer>
        void construct_samples(const RandomAccessContainer& sa, const unsigned char* str);

        void copy(const csa_sada& csa) {
            m_psi = csa.m_psi;
            m_sa_sample = csa.m_sa_sample;
            m_isa_sample = csa.m_isa_sample;
            m_char2comp  = csa.m_char2comp;
            m_comp2char  = csa.m_comp2char;
            m_C = csa.m_C;
            m_sigma		 = csa.m_sigma;
            m_psi_wrapper = psi_type(this);
            m_bwt = bwt_type(this);
        };

    public:
        const int_vector<8>& char2comp;
        const int_vector<8>& comp2char;
        const int_vector<64>& C;
        const uint16_t& sigma;
        const psi_type& psi;
        const bwt_type& bwt;
        const sa_sample_type& sa_sample;
        const isa_sample_type& isa_sample;


        //! Default Constructor
        csa_sada():char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma) ,psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample) {
            m_psi_buf = NULL;
            m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];
            util::assign(m_char2comp, int_vector<8>(256, 0));
            util::assign(m_comp2char, int_vector<8>(256, 0));
            util::assign(m_C, int_vector<64>(257, 0));
            m_sigma = 0;
        }
        //! Default Destructor
        ~csa_sada() {
            if (m_psi_buf != NULL)
                delete [] m_psi_buf;
        }
        //! Copy constructor
        csa_sada(const csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>& csa):char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma), psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample) {
            m_psi_buf = NULL;
            m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];
            copy(csa);
        }

        //! Construct the csa_sada from another (compressed) suffix array and the original text
        /*! \param sa  The existing (compressed) suffix array for the string str.
         *	\param str The text for which sa was build.
         */
        template<typename RandomAccessContainer>
        csa_sada(const RandomAccessContainer& sa, const unsigned char* str);

        //! Construct the csa_sada from another (compressed) suffix array and the original text
        /*! \param sa  The existing (compressed) suffix array for the string str.
         *	\param str The text for which the suffix array sa was build.
         */
        template<typename RandomAccessContainer>
        csa_sada(RandomAccessContainer& sa, const unsigned char* str);

        //! Construct the csa_sada from the int_vector_file_buffers of the text, the suffix array and the inverse suffix array
        template<class size_type_class, uint8_t int_width, class size_type_class_1>
        csa_sada(int_vector_file_buffer<8, size_type_class>& bwt_buf,
                 int_vector_file_buffer<int_width, size_type_class_1>& sa_buf);

        csa_sada(tMSS& file_map, const std::string& dir, const std::string& id);

        void construct(tMSS& file_map, const std::string& dir, const std::string& id);

        //! Constructor for the CSA taking a string for that the CSA should be calculated
        /*!	\param str The text for which the CSA should be constructed.
         */
        csa_sada(const unsigned char* str);

        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const {
            return m_psi.size();
        }

        //! Returns the largest size that csa_sada can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return EncVector::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_psi.empty();
        }

        //! Swap method for csa_sada
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param csa csa_sada to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>& csa);

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const;

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const;

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Required for the STL Random Access Container Concept.
         * \par Time complexity
         *      \f$ \Order{s_{SA}\cdot t_{\Psi}} \f$, where every \f$s_{SA}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline value_type operator[](size_type i)const;

        //! ()-operator return inverse suffix array values
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         *   \par Time complexity
         *      \f$ \Order{s_{SA^{-1}}\cdot t_{\Psi}} \f$, where every \f$s_{SA^{-1}}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline value_type operator()(size_type i)const;

        //! Assignment Operator.
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        csa_sada& operator=(const csa_sada& csa);

        //! Equality Operator
        /*! Two Instances of csa_sada are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const csa_sada& csa)const;

        //! Unequality Operator
        /*! Two Instances of csa_sada are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const csa_sada& csa)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

        uint32_t get_sample_dens() const;

        uint32_t get_psi_sample_dens() const;
        void set_psi_sample_dens(const uint32_t sample_dens);

        //! Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\returns The number of occurences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *		\f$ \Order{\log n t_{\Psi}} \f$
         */
        size_type rank_bwt(size_type i, const unsigned char c)const {
            unsigned char cc = m_char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return 0
                return 0;
            if (i == 0)
                return 0;
            assert(i <= size());

            size_type lower_b, upper_b; // lower_b inclusive, upper_b exclusive

            const size_type sd = m_psi.get_sample_dens();
            size_type lower_sb = (m_C[cc]+sd-1)/sd; // lower_sb inclusive
            size_type upper_sb = (m_C[cc+1]+sd-1)/sd; // upper_sb exclusive
            while (lower_sb+1 < upper_sb) {
                size_type mid = (lower_sb+upper_sb)/2;
                if (m_psi.sample(mid) >= i)
                    upper_sb = mid;
                else
                    lower_sb = mid;
            }

            if (lower_sb == upper_sb) { // the interval was smaller than sd
                lower_b = m_C[cc]; upper_b = m_C[cc+1];
            } else if (lower_sb > (m_C[cc]+sd-1)/sd) { // main case
                lower_b = lower_sb*sd;
                uint64_t* p=m_psi_buf;
                // extract the psi values between two samples
                m_psi.get_inter_sampled_values(lower_sb, p);
                p = m_psi_buf;
                uint64_t smpl = m_psi.sample(lower_sb);
                // handle border cases
                if (lower_b + m_psi.get_sample_dens() >= m_C[cc+1])
                    m_psi_buf[ m_C[cc+1]-lower_b ] = size()-smpl;
                else
                    m_psi_buf[ m_psi.get_sample_dens() ] = size()-smpl;
                // search the result linear
                while ((*p++)+smpl < i);

                return p-1-m_psi_buf + lower_b - m_C[cc];
            } else { // lower_b == (m_C[cc]+sd-1)/sd and lower_sb < upper_sb
                if (m_psi.sample(lower_sb) >= i) {
                    lower_b = m_C[cc];
                    upper_b = lower_sb * sd + 1;
                } else {
                    lower_b = lower_sb * sd;
                    upper_b = std::min(upper_sb*sd, m_C[cc+1]);
                }
            }

            // binary search the interval [C[cc]..C[cc+1]-1] for the result
//			size_type lower_b = m_C[cc], upper_b = m_C[cc+1]; // lower_b inclusive, upper_b exclusive
            while (lower_b+1 < upper_b) {
                size_type mid = (lower_b+upper_b)/2;
                if (m_psi[mid] >= i)
                    upper_b = mid;
                else
                    lower_b = mid;
            }
            if (lower_b > m_C[cc])
                return lower_b - m_C[cc] + 1;
            else { // lower_b == m_C[cc]
                return m_psi[lower_b] < i;// 1 if m_psi[lower_b]<i, 0 otherwise
            }
        }

        //! Calculates the ith occurence of symbol c in the BWT of the original text.
        /*!
         *"  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *	\returns The ith occurence symbol c in the the BWT or size() if no ith occurence of the symbol exists in the BWT.
         *  \par Time complexity
         *		\f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const unsigned char c)const {
            assert(i > 0);
            unsigned char cc = m_char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return 0
                return size();
            assert(cc != 255);
            if (m_C[cc]+i-1 <  m_C[cc+1]) {
                return m_psi[m_C[cc]+i-1];
            } else
                return size();
        }

#ifdef MEM_INFO
        //! Print some infos about the size of the compressed suffix tree
        void mem_info(std::string label="")const {
            if (label=="")
                label = "csa";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<"\n,";
            m_psi.mem_info("$\\Psi$");
            std::cout<<",";
            sa_sample.mem_info("sa_sample");
            std::cout<<",";
            isa_sample.mem_info("isa_sample");
            std::cout << ")\n";
        }
#endif

};

// == template functions ==

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::csa_sada(const unsigned char* str):char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma), psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample)
{
    m_psi_buf = NULL;
    m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];
    csa_uncompressed sa(str);
//	size_type n = strlen((const char*)str);
//	int_vector<> sa(n+1, 0, bit_magic::l1BP(n+1)+1);
//	algorithm::calculate_sa(str, n+1, sa);	 // calculate the suffix array sa of str
//	assert(sa.size() == n+1);
    algorithm::set_text<csa_sada>(str, sa.size(), m_C, m_char2comp, m_comp2char, m_sigma);
    construct(sa, str);
//	if( n+1 > 0 and n+1 != size() )
//		throw std::logic_error("csa_sada: text size differ with sa size!");
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
template<typename RandomAccessContainer>
csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::csa_sada(const RandomAccessContainer& sa, const unsigned char* str):char2comp(m_char2comp), comp2char(m_comp2char),C(m_C), sigma(m_sigma), psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample)
{
    m_psi_buf = NULL;
    m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];

    size_type n = 1;
    if (str != NULL) {
        n = strlen((const char*)str);
    }
    algorithm::set_text<csa_sada>(str, n+1, m_C, m_char2comp, m_comp2char, m_sigma);
    assert(sa.size() == n+1);
    construct(sa, str);
    if (n+1 > 0 and n+1 != size())
        throw std::logic_error(util::demangle(typeid(this).name())+": text size differ with sa size!");
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
template<typename RandomAccessContainer>
csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::csa_sada(RandomAccessContainer& sa, const unsigned char* str):char2comp(m_char2comp), comp2char(m_comp2char),C(m_C), sigma(m_sigma), psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample)
{
    m_psi_buf = NULL;
    m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];

    size_type n = 1;
    if (str != NULL) {
        n = strlen((const char*)str);
    }
    assert(sa.size() == n+1);
    algorithm::set_text<csa_sada>(str, n+1, m_C, m_char2comp, m_comp2char, m_sigma);
    construct(sa, str);
    if (n+1 > 0 and n+1 != size())
        throw std::logic_error("csa_sada: text size differ with sa size!");
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
template<class size_type_class, uint8_t int_width, class size_type_class_1>
csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::csa_sada(int_vector_file_buffer<8, size_type_class>& bwt_buf,
        int_vector_file_buffer<int_width, size_type_class_1>& sa_buf):
    char2comp(m_char2comp), comp2char(m_comp2char),C(m_C), sigma(m_sigma), psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample)
{
    m_psi_buf = NULL;
    m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];
    bwt_buf.reset(); sa_buf.reset();
    size_type n = bwt_buf.int_vector_size;
    algorithm::set_text<csa_sada>(bwt_buf, n, m_C, m_char2comp, m_comp2char, m_sigma);
    assert(sa_buf.int_vetor_size == sa_buf.int_vector_size);

    size_type cnt_chr[256] = {0};
    for (uint32_t i=0; i<m_sigma; ++i)
        cnt_chr[m_comp2char[i]] = C[i];
    stop_watch sw; sw.start();
//	write_R_output(sw, "psi function", "construct","begin");
    // calculate psi
    {
        bwt_buf.reset();
        int_vector<> temp(n, 0, bit_magic::l1BP(n)+1);
        for (size_type i=0, r_sum=0, r=bwt_buf.load_next_block(); r_sum < n;) {
            for (; i < r_sum+r; ++i) {
                temp[ cnt_chr[ bwt_buf[i-r_sum] ]++ ] = i;
            }
            r_sum += r; r = bwt_buf.load_next_block();
        }
        util::store_to_file(temp, "/tmp/deleteme");
    }
//	write_R_output(sw, "psi function", "construct","end");
    int_vector_file_buffer<> psi_buf("/tmp/deleteme");
//	write_R_output(sw, "encoded psi", "construct", "begin");
    m_psi = EncVector(psi_buf);
//	write_R_output(sw, "encoded psi", "construct", "end");
    std::remove("/tmp/deleteme");
    algorithm::set_sa_and_isa_samples<csa_sada>(sa_buf, m_sa_sample, m_isa_sample);
    m_psi_wrapper = psi_type(this);
    m_bwt = bwt_type(this);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::csa_sada(tMSS& file_map, const std::string& dir, const std::string& id):
    char2comp(m_char2comp), comp2char(m_comp2char),C(m_C), sigma(m_sigma), psi(m_psi_wrapper), bwt(m_bwt), sa_sample(m_sa_sample), isa_sample(m_isa_sample)
{
    construct(file_map, dir, id);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::construct(tMSS& file_map, const std::string& dir, const std::string& id)
{
    if (m_psi_buf = NULL) {
        m_psi_buf = new uint64_t[m_psi.get_sample_dens()+1];
    }
    if (file_map.find("bwt") == file_map.end()) { // if bwt is not already stored on disk => construct bwt
        construct_bwt(file_map, dir, id);
    }
    int_vector_file_buffer<8> bwt_buf(file_map["bwt"].c_str());
    size_type n = bwt_buf.int_vector_size;
    algorithm::set_text<csa_sada>(bwt_buf, n, m_C, m_char2comp, m_comp2char, m_sigma);

    size_type cnt_chr[256] = {0};
    for (uint32_t i=0; i<m_sigma; ++i)
        cnt_chr[m_comp2char[i]] = C[i];
    stop_watch sw; sw.start();
    write_R_output("csa", "construct PSI","begin",1,0);
    // calculate psi
    {
        bwt_buf.reset();
        int_vector<> psi(n, 0, bit_magic::l1BP(n)+1);
        for (size_type i=0, r_sum=0, r=bwt_buf.load_next_block(); r_sum < n;) {
            for (; i < r_sum+r; ++i) {
                psi[ cnt_chr[ bwt_buf[i-r_sum] ]++ ] = i;
            }
            r_sum += r; r = bwt_buf.load_next_block();
        }
        if (!util::store_to_file(psi, (dir+"psi_"+id).c_str())) {
            throw std::ios_base::failure("#csa_sada: Cannot store PSI to file system!");
        } else {
            file_map["psi"] = dir+"psi_"+id;
        }
    }
    write_R_output("csa", "construct PSI","end",1,0);
    int_vector_file_buffer<> psi_buf(file_map["psi"].c_str());
//	write_R_output(sw, "encoded psi", "construct", "begin");
    m_psi = EncVector(psi_buf);
//	write_R_output(sw, "encoded psi", "construct", "end");
    int_vector_file_buffer<>  sa_buf(file_map["sa"].c_str());
    algorithm::set_sa_and_isa_samples<csa_sada>(sa_buf, m_sa_sample, m_isa_sample);
    m_psi_wrapper = psi_type(this);
    m_bwt = bwt_type(this);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
uint32_t csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::get_sample_dens()const
{
    return SampleDens;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
uint32_t csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::get_psi_sample_dens()const
{
    return m_psi.get_sample_dens();
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::set_psi_sample_dens(const uint32_t sample_dens)
{
    m_psi.get_sample_dens();
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
template<typename RandomAccessContainer>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::construct_samples(const RandomAccessContainer& sa, const unsigned char* str)
{
    m_sa_sample.set_int_width(bit_magic::l1BP(sa.size())+1);
    m_sa_sample.resize((sa.size()+get_sample_dens()-1)/get_sample_dens());
    typename RandomAccessContainer::size_type i=0, idx=0;
    for (typename RandomAccessContainer::const_iterator it = sa.begin(); i < sa.size(); it += (difference_type)get_sample_dens(), i += get_sample_dens(), ++idx) {
        m_sa_sample[idx] = *it;
    }
//	const uint32_t SampleDens, uint32_t InvSampleDensISA = get_sample_dens()*16;
    m_isa_sample.set_int_width(bit_magic::l1BP(sa.size())+1);
    m_isa_sample.resize((sa.size()+(InvSampleDens)-1)/(InvSampleDens));
    i = 0;
    for (typename RandomAccessContainer::const_iterator it = sa.begin(), end = sa.end(); it != end; ++it, ++i) {
        if ((*it % InvSampleDens) == 0) {
            m_isa_sample[*it/InvSampleDens] = i;
        }
    }
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
template<typename RandomAccessContainer>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::construct(const RandomAccessContainer& sa, const unsigned char* str)
{
    construct_samples(sa, str);
#ifdef SDSL_DEBUG
    std::cerr<<"create encoded psi"<<std::endl;
#endif
    if (sa.psi.size() > 0) {
        m_psi = EncVector(sa.psi);
    } else {
        m_psi = EncVector();
    }
    m_psi_wrapper = psi_type(this);
    m_bwt = bwt_type(this);
#ifdef SDSL_DEBUG
    std::cerr<<"encoded psi created"<<std::endl;
#endif
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
template<typename RandomAccessContainer>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::construct(RandomAccessContainer& sa, const unsigned char* str)
{
    construct_samples(sa, str);
#ifdef SDSL_DEBUG
    std::cerr<<"create encoded psi"<<std::endl;
#endif
    if (sa.psi.size() > 0) {
        m_psi = EncVector(sa.psi);
    } else {
        m_psi = EncVector();
    }
    m_psi_wrapper = psi_type(this);
    m_bwt = bwt_type(this);
    {
        // clear sa
        RandomAccessContainer tmp;
        tmp.swap(sa);
    }
#ifdef SDSL_DEBUG
    std::cerr<<"encoded psi and created; sa destroyed"<<std::endl;
#endif
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
typename csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::const_iterator csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::begin()const
{
    return const_iterator(this, 0);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
typename csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::const_iterator csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::end()const
{
    return const_iterator(this, size());
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
inline typename csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::value_type csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::operator[](size_type i)const
{
    size_type off = 0;
    while (i % SampleDens) {// while i mod SampleDens != 0 (SA[i] is not sampled)   SG: auf keinen Fall get_sample_dens nehmen, ist total langsam
        i = m_psi[i];       // go to the position where SA[i]+1 is located
        ++off;              // add 1 to the offset
    }
    value_type result = m_sa_sample[i/SampleDens];
    if (result < off) {
        return m_psi.size()-(off-result);
    } else
        return result-off;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
inline typename csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::value_type csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::operator()(size_type i)const
{
    value_type result = m_isa_sample[i/InvSampleDens]; // get the rightmost sampled isa value
    i = i % InvSampleDens;
    while (i--) {
        result = m_psi[result];
    }
//	assert(((*this)[result])==j);
    return result;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
csa_sada<EncVector,SampleDens, InvSampleDens, fixedIntWidth>& csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::operator=(const csa_sada<EncVector,SampleDens, InvSampleDens, fixedIntWidth>& csa)
{
    if (this != &csa) {
        copy(csa);
    }
    return *this;
}


template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
typename csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::serialize(std::ostream& out)const
{
    size_type written_bytes = 0;
    written_bytes += m_psi.serialize(out);
    written_bytes += m_sa_sample.serialize(out);
    written_bytes += m_isa_sample.serialize(out);
    written_bytes += m_char2comp.serialize(out);
    written_bytes += m_comp2char.serialize(out);
    written_bytes += m_C.serialize(out);
    written_bytes += util::write_member(m_sigma, out);
    return written_bytes;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::load(std::istream& in)
{
    m_psi.load(in);
    m_sa_sample.load(in);
    m_isa_sample.load(in);
    m_char2comp.load(in);
    m_comp2char.load(in);
    m_C.load(in);
    util::read_member(m_sigma, in);
    m_psi_wrapper = psi_type(this);
    m_bwt = bwt_type(this);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
void csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::swap(csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>& csa)
{
    if (this != &csa) {
        m_psi.swap(csa.m_psi);
        m_sa_sample.swap(csa.m_sa_sample);
        m_isa_sample.swap(csa.m_isa_sample);
        m_char2comp.swap(csa.m_char2comp);
        m_comp2char.swap(csa.m_comp2char);
        m_C.swap(csa.m_C);
        std::swap(m_sigma, csa.m_sigma);
        m_psi_wrapper = psi_type(this);
        csa.m_psi_wrapper = psi_type(&csa);
        m_bwt = bwt_type(this);
        csa.m_bwt = bwt_type(&csa);
    }
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
bool csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::operator==(const csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>& csa)const
{
    for (uint16_t i=0; i<256; ++i)
        if (m_char2comp[i] != csa.m_char2comp[i] or m_comp2char[i] != csa.m_comp2char[i] or m_C[i] != csa.m_C[i])
            return false;
    return m_psi == csa.m_psi and m_sa_sample == csa.m_sa_sample and m_isa_sample == csa.m_isa_sample and m_C[256] == csa.m_C[256] and m_sigma == csa.m_sigma;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
bool csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::operator!=(const csa_sada<EncVector, SampleDens, InvSampleDens, fixedIntWidth>& csa)const
{
    return !(*this == csa);
//	return m_psi != csa.m_psi or m_sa_sample != csa.m_sa_sample or m_isa_sample != csa.m_isa_sample;
}

} // end namespace sdsl

#endif
