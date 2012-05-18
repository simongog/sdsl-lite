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
/*! \file enc_vector_dna.hpp
   \brief enc_vector_dna.hpp contains the sdsl::enc_vector_dna class.
   \author Simon Gog
*/
#include "int_vector.hpp"
#include "fibonacci_coder.hpp"
#include "iterators.hpp"

#ifndef SDSL_ENC_VECTOR_DNA
#define SDSL_ENC_VECTOR_DNA

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t fixedIntWidth>
struct enc_vector_dna_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct enc_vector_dna_trait<32> {
    typedef int_vector<32> int_vector_type;
};

//! An immutable space-saving vector class for unsigned positiv integers of the psi-values of dna data.
/*! It encodes each integer with its fibonacci-code and still provides constant time access.
 *  The values of a enc_vector_dna are immutable after the constructor call. The class
 *  could be parametrized with sample denisty.
 *
 * @ingroup int_vector
 */
template<uint32_t SampleDens = 8, uint8_t fixedIntWidth=0>
class enc_vector_dna
{
    public:
        typedef uint64_t 							value_type;  	// STL Container requirement
        typedef random_access_const_iterator<enc_vector_dna> iterator;// STL Container requirement
        typedef iterator							const_iterator; // STL Container requirement
        typedef const value_type		 			reference;
        typedef const value_type 					const_reference;
//	typedef SDSBitVectorReference*				pointer;
        typedef const value_type*					const_pointer;
        typedef ptrdiff_t 							difference_type;// STL Container requirement
        typedef int_vector<>::size_type				size_type;		// STL Container requirement
        typedef coder::fibonacci					coder;
        static  const uint32_t 						sample_dens	= SampleDens;

        int_vector<0> 	m_z; 		// compressed bit stream
    private:
        typename enc_vector_dna_trait<fixedIntWidth>::int_vector_type   m_sample_vals_and_pointer;
        size_type		m_elements;    // number of elements

        // workaround function for the constructor
        void construct() {
            m_elements = 0;
        }
        void copy(const enc_vector_dna& v);

        void write_sumblock(uint64_t x, uint64_t blocknr, uint64_t*& sumblock, bool*& overflow) {
            if (!overflow[blocknr]) {
                if (x > 255) {
                    overflow[blocknr] = true; sumblock[blocknr] = 0;
                } else {
                    sumblock[blocknr] += x;
                    if (sumblock[blocknr]>255) {
                        overflow[blocknr] = true;
                        sumblock[blocknr] = 0;
                    }
                }
            }

        }

    public:
        //! Default Constuctor
        enc_vector_dna() {
            construct();
        };
        //! Copy constructor
        /*! \param v The enc_vector_dna to copy.
          	Required for the Assignable Concept of the STL
         */
        enc_vector_dna(const enc_vector_dna& v);

        //! Constructor for a Container of positiv integers.
        /*! \param c A container of positive integers.
            \par The container is used to build the EncVector of the
        	     integer sequence.
          */
        template<class Container>
        enc_vector_dna(const Container& c) {
            construct();
            init(c);
        };

        template<class Container>
        void init(const Container& c);

        //! Default Destructor
        ~enc_vector_dna() {
        };

        //! The number of elements in the enc_vector_dna.
        /*!
         	Required for the Container Concept of the STL.
        	\sa max_size
         */
        size_type size()const;

        //! Return the largest size that this container can ever have.
        /*! Required for the Container Concept of the STL.
         */
        static size_type max_size();

        //!	Returns if the enc_vector_dna is empty.
        /*! Equivalent to size() == 0.
         *
         * 	Required for the STL Container Concept.
         *  \sa size()
         */
        bool empty() const;

        //! Swap method for enc_vector_dna
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param v enc_vector_dna to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(enc_vector_dna& v);

        //! Iterator that points to the first element of the enc_vector_dna.
        /*!
         * 	Required for the Container Concept of the STL.
         *  \sa end()
         */
        const const_iterator begin()const;

        //! Iterator that points to the position after the last element of the enc_vector_dna.
        /*!
         *	Required for the Container Concept of the STL
         *  \sa begin()
         */
        const const_iterator end()const;

        // Iterator that points to the last element of the enc_vector_dna.
        /*
         * 	Required for the Container Concept of the STL.
         *  \sa rend()
         */
//		reverse_iterator rbegin()const;

        // Iterator that points to the position before the first element of the enc_vector_dna.
        /*
         *	Required for the Container Concept of the STL
         *  \sa rbegin()
         */
//		reverse_iterator rend()const;

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         *
         *  Required for the STL Random Access Container Concept.
         */
        value_type operator[](size_type i)const;

        //! Assignment Operator
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        enc_vector_dna& operator=(const enc_vector_dna& v);

        //! Equality Operator
        /*! Two enc_vector_dnas are equal if all member variables are equal
         *  (including the sample density of the enc_vector_dnas).
         *  \note If the sample density is not equal you should use
         *  SDSAlgorithm::equal_container_values to compare two enc_vector_dnas.
         *
         * 	Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const enc_vector_dna& v)const;

        //! Unequality Operator
        /*! Two enc_vector_dnas are unuequal if not all member variables are equal
         *  (including the sample density of the enc_vector_dnas).
         *  \note If the sample density is not equal you should use
         *  SDSAlgorithm::equal_container_values to compare two enc_vector_dnas.
         *
         * 	Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const enc_vector_dna& v)const;

        //! Serialzes the enc_vector_dna to a stream.
        /*! \param out Outstream to write the data structure.
            \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load the enc_vector_dna from a stream.
        void load(std::istream& in);

        //! Returns the ith sample of enc_vector_dna
        /*! \param i The index of the sample. 0 <= i < size()/SampleDens
         *  \return The value of the ith sample.
         */
        value_type sample(const size_type i) const;
};


template<uint32_t SampleDens, uint8_t fixedIntWidth>
inline typename enc_vector_dna<SampleDens,fixedIntWidth>::value_type enc_vector_dna<SampleDens,fixedIntWidth>::sample(const size_type i)const
{
    if (i+1 == 0 || i >= m_elements/SampleDens) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_dna::sample[](size_type); i >= size()!");
        return 0;
    }
    return m_sample_vals_and_pointer[i<<1];
}


inline uint64_t decode64bit(uint64_t w)
{
    uint64_t result = 0;
    while (w) {
        uint16_t temp = coder::fibonacci::Fib2bin_0_16_greedy[w&0xFFFF], shift;
        if ((shift=(temp>>11)) > 0) {
            result += (temp & 0x7FFULL);
            w >>= shift;
        } else {
            temp = 0;
            do {
                result += coder::fibonacci::Fib2bin_0_95[(temp<<12) | (w&0xFFF)];
                if ((shift =  coder::fibonacci::Fib2binShift[w&0x1FFF]) > 0) {
                    w >>= shift;
                    break;
                } else {
                    w >>= 12;
                    temp++;
                }
            } while (1);
        }
    }
    return result;
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
inline typename enc_vector_dna<SampleDens,fixedIntWidth>::value_type enc_vector_dna<SampleDens,fixedIntWidth>::operator[](size_type i)const
{
    if (i+1 == 0 || i >= m_elements) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_dna::operator[](size_type); idx >= size()!");
        return 0;
    }
    const size_type idx	 = i/SampleDens;
    i 				-= SampleDens*idx; // values to decode
    uint64_t result = m_sample_vals_and_pointer[idx<<1];
    if (!i) // if i==0
        return result;
//	uint64_t expected = result + coder::fibonacci::decode1<true, false, int*>(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], i);

    uint8_t offset			= m_sample_vals_and_pointer[(idx<<1)+1] & 0x3F;
    const uint64_t* data	= m_z.data() + (m_sample_vals_and_pointer[(idx<<1)+1] >> 6);
    uint64_t w				= (*data) & ~bit_magic::Li1Mask[offset];
//	std::cerr<<"0 w="<<w<<std::endl;
    bool carry = false;
    uint64_t m = bit_magic::all11BPs(w, carry); // calc end positions of fibonacci encoded numbers
    size_type cnt = bit_magic::b1Cnt(m); // count ends of fibonacci encoded numbers

    if (cnt >= i) {
        // decode i values in a block and return result
        w &= bit_magic::Li1Mask[bit_magic::i1BP(m,i)+1];
//		assert( result + decode64bit(w>>offset) == expected );
        return result + decode64bit(w>>offset);
    } else {
        // 0 < cnt < i
//		if(cnt==0) goto slow_decode;
        uint8_t last1pos = bit_magic::l1BP(m)+/*1*/(cnt>0),temp=0;
        w &= bit_magic::Li1Mask[last1pos];
        result += decode64bit(w>>offset);
//	uint64_t expected1 = m_sample_vals_and_pointer[idx<<1] + coder::fibonacci::decode1<true, false, int*>(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], cnt);
//		assert(result == expected1);

        uint8_t  blocknr	= (data-m_z.data())%9;
        uint64_t sumblock	= (*(data + 8 - blocknr))>>(blocknr<<3);
        do {
            ++data; ++blocknr;
            sumblock >>= 8;
            if (blocknr==8) {
                ++data;
                sumblock=*(data+8);
                blocknr=0;
            }
            m = bit_magic::all11BPs(*data,carry);
            cnt += (temp=bit_magic::b1Cnt(m));
            if (sumblock&0xFF) {
                result += sumblock&0xFF;
//				if(cnt<=i){
//	uint64_t expected2 = m_sample_vals_and_pointer[idx<<1] + coder::fibonacci::decode1<true, false, int*>(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], cnt);
//				assert(result==expected2);
//				}
            } else
                goto slow_decode;
        } while (cnt < i);
        if (cnt==i) {
//			assert(result == expected);
            return result;
        }
        offset = bit_magic::i1BP(m,temp+i-cnt)+1;
        w = (*data) >> offset;
        w &= bit_magic::Li1Mask[bit_magic::l1BP(m>>offset)+1];
//		assert( result == expected );
        return result-decode64bit(w);
    }

slow_decode:

    result = m_sample_vals_and_pointer[idx<<1];
    return result + coder::fibonacci::decode1<true, false, int*>(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], i);
}
/*
template<uint32_t SampleDens, uint8_t fixedIntWidth>
inline typename enc_vector_dna<SampleDens,fixedIntWidth>::value_type enc_vector_dna<SampleDens,fixedIntWidth>::operator[](const size_type i)const {
	if( i+1 == 0 || i >= m_elements  ){
		throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_dna::operator[](size_type); idx >= size()!");
		return 0;
	}
	size_type idx = i/SampleDens;
	size_type n   = i-SampleDens*idx; // values to decode
	uint64_t result = m_sample_vals_and_pointer[idx<<1];
	if(n==0)
		return result;
//if(n==69)
//std::cout<<"n="<<n<<std::endl;
	// n > 0
	int16_t offset = m_sample_vals_and_pointer[(idx<<1)+1] & 0x3F;
	uint64_t w = m_sample_vals_and_pointer[(idx<<1)+1] >> 6;
	const uint64_t *data = m_z.data() + w;
	uint8_t blocknr = w%9;
//std::cout<<"blocknr="<<(int)blocknr<<std::endl;

//uint64_t checkbuf[1024] = {0};
//checkbuf[0]= *data;
//uint64_t checkres = m_sample_vals_and_pointer[idx<<1];
//uint64_t *cb=checkbuf;
//int16_t checkoff = offset;

	w = *data;
	bool carry = false;
	uint64_t m = bit_magic::all11BPs(w&~bit_magic::Li1Mask[offset], carry); // calc end positions of fibonacci encoded numbers
	size_type cnt = bit_magic::b1Cnt(m); // count ends of fibonacci encoded numbers
	if( cnt >= n  ){
//std::cout<<"Case a: result so far = "<<result<<" cnt="<<cnt<<" n="<<n<<" offset="<<(int)offset<<std::endl;
//std::cout<<"w="<<w<<std::endl;
//std::cout<<"m="<<m<<std::endl;
		// decode n values in a block and return result
		w >>= offset;
		m >>= offset;
		uint16_t to_decode = bit_magic::i1BP(m,n)+1;
//std::cout<<"to_decode = "<<to_decode<<std::endl;
		w &= bit_magic::Li1Mask[to_decode];
//std::cout<<"w="<<w<<std::endl;
		// decode n values in a block and return result
//std::cout<<"decoding "<<decode64bit(w)<<std::endl;
		return result + decode64bit(w);
	}else{
//std::cout<<"offset="<<(int)offset<<std::endl;
//std::cout<<"cnt="<<cnt<<std::endl;
//std::cout<<"w="<<w<<std::endl;
//std::cout<<"m="<<m<<" w="<<w<<" cnt="<<cnt<<" offset="<<offset<<std::endl;
		m >>= offset;
		w >>= offset;
		if(m){
			result += decode64bit( w & bit_magic::Li1Mask[bit_magic::l1BP(m)+1] );
		}
		m <<= offset;
		if(offset) m |= (1ULL<<(offset-1));
		w <<= offset;
		uint64_t sumblock = *(data + 8 - blocknr), oldm;
		sumblock >>= (blocknr<<3);
		uint16_t temp;
		do{
			++data; ++blocknr;
			sumblock >>= 8;
//std::cout<<"sumblock="<<sumblock<<std::endl;
			if(blocknr==8){sumblock=*data; ++data; blocknr=0;}
// *(++cb)=*data;
			oldm = m;
			m    = bit_magic::all11BPs(*data, carry);
			cnt += (temp=bit_magic::b1Cnt(m));
			if( cnt < n ){
				if( (sumblock & 0xFF) ){// sum is precalculated
//if(n==69)
//std::cout<<"Case 1: sum = "<< (sumblock&0xFF) <<" result so far = "<<result<<" cnt="<<cnt<<" n="<<n<<" w="<<*data<<std::endl;
					result += (sumblock & 0xFF);
				}
				else{// sum is not precalculated
					offset = bit_magic::l1BP(oldm)+1; // find the leftmost position where an encoded value terminates in the previous word
//if(n==69)
//std::cout<<"Case 2: sum not precalculated "<<" cnt="<<cnt<<" n="<<n<<" offset="<<(int)offset<<std::endl;
					if(m == 0){// a big value ends in the next word => blocksum of next word = 0, and cnt < n
//std::cout<<"Case 2a: m==0 "<<std::endl;
						uint64_t buffer[3];
						buffer[0] = *(data-1-(blocknr==0));  // get the last word
						buffer[1] = *data;
//						oldcarry = carry;
						++data; ++blocknr;
						sumblock >>= 8;
						if(blocknr==8){sumblock=*data; ++data; blocknr=0;}
// *(++cb)=*data;
						buffer[2] = *data;
						oldm = m;
						m    = bit_magic::all11BPs(*data, carry);
//std::cout<<"Carry = "<<(int)carry<<std::endl;
//std::cout<<buffer[0]<<std::endl;
//std::cout<<buffer[1]<<std::endl;
//std::cout<<buffer[2]<<std::endl;
//std::cout<<"m="<<m<<std::endl;
						cnt += (temp=bit_magic::b1Cnt(m));
//std::cout<<"temp="<<(int)temp<<std::endl;
						if(cnt >= n){
//std::cout<<"decode = "<< coder::fibonacci::decode<true,false, int*>(buffer, offset, n + temp - cnt) << std::endl << "result = " << result << std::endl;
//std::cout<<"decode check = "<< coder::fibonacci::decode<true,false, int*>(checkbuf, checkoff, n) << std::endl;
//std::cout<<"decode check = "<< coder::fibonacci::decode<true,false, int*>(checkbuf, checkoff, n-1) << std::endl;
							return result + coder::fibonacci::decode<true, false, int*>(buffer, offset, n + temp - cnt);
						}else{
							result += coder::fibonacci::decode<true, false, int*>(buffer, offset, temp);
						}
					}else{
//if(n==69)
//std::cout<<"Case 2b: "<<std::endl;
						w = *data;
						uint16_t to_decode = bit_magic::l1BP(m)+1;
						to_decode += (64-offset);
//						std::cout<<"to_decode"<<(int)to_decode<<std::endl;
						if( to_decode <= 64 ){
//if(n==69)
//std::cout<<"Case 2b1: "<<std::endl;
//							std::cout<<"w="<<w<<std::endl;
//							std::cout<<"m="<<m<<std::endl;
							w <<= (64-offset);
//							std::cout<<"w="<<w<<std::endl;
							w |= ((*(data-1-(blocknr==0)) >> offset) & bit_magic::Li1Mask[64-offset]);
//							std::cout<<"last w="<<(*(data-1-(blocknr==0)))<<std::endl;
							w &= bit_magic::Li1Mask[to_decode];
//							std::cout<<"w="<<w<<std::endl;
							result += decode64bit(w);
						}else{
//if(n==69)
//std::cout<<"Case 2b2: "<<std::endl;
							uint64_t buffer[2] = { *(data-1-(blocknr==0)), *data  };// falsch wenn buffer[0] das ganz erste wort ist
							result += coder::fibonacci::decode<true, false, int*>(buffer, offset, temp);
						}
					}
				}
			}else{// cnt >= n, => w contains at least one 1 bit
//if(n==69)
//std::cout<<"Case 3:"<<std::endl;
				uint16_t to_decode = bit_magic::i1BP(m, temp)+1;
				w = *data & bit_magic::Li1Mask[to_decode];
				offset = bit_magic::l1BP(oldm)+1;
				to_decode += (64-offset);
				if( to_decode <= 64 ){
					w <<= (64-offset); // space for the encoded value from the previous word
					w |= ((*(data-1-(blocknr==0)) >> offset) & bit_magic::Li1Mask[64-offset]);
					return result + decode64bit(w);
				}else{
					uint64_t buffer[2] = { *(data-1-(blocknr==0)), *data  };
					return result + coder::fibonacci::decode<true, false, int*>(buffer, offset, temp);
				}
			}
		}while(1);
	}
	return result;
	//return m_sample_vals_and_pointer[idx<<1] + coder::decode_prefix_sum(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], i-SampleDens*idx );
}
*/

template<uint32_t SampleDens, uint8_t fixedIntWidth>
inline enc_vector_dna<>::size_type enc_vector_dna<SampleDens,fixedIntWidth>::size()const
{
    return m_elements;
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
inline enc_vector_dna<>::size_type enc_vector_dna<SampleDens,fixedIntWidth>::max_size()
{
    return int_vector<>::max_size()/2; // each element could possible occupy double space with selfdelimiting codes
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
inline bool enc_vector_dna<SampleDens,fixedIntWidth>::empty()const
{
    return 0==m_elements;
}


template<uint32_t SampleDens, uint8_t fixedIntWidth>
void enc_vector_dna<SampleDens,fixedIntWidth>::copy(const enc_vector_dna<SampleDens,fixedIntWidth>& v)
{
    m_z					= v.m_z;				// copy compressed bit stream
    m_sample_vals_and_pointer		= v.m_sample_vals_and_pointer;      // copy sample values
    m_elements			= v.m_elements;			// copy number of stored elements
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
enc_vector_dna<SampleDens,fixedIntWidth>::enc_vector_dna(const enc_vector_dna& v)
{
    copy(v);
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
enc_vector_dna<SampleDens,fixedIntWidth>& enc_vector_dna<SampleDens,fixedIntWidth>::operator=(const enc_vector_dna<SampleDens,fixedIntWidth>& v)
{
    if (this != &v) {// if v and _this_ are not the same object
        copy(v);
    }
    return *this;
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
bool enc_vector_dna<SampleDens,fixedIntWidth>::operator==(const enc_vector_dna<SampleDens,fixedIntWidth>& v)const
{
    if (this == &v)
        return true;
    return	 	m_elements == v.m_elements
                and	m_z == v.m_z
                and	m_sample_vals_and_pointer == v.m_sample_vals_and_pointer;
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
bool enc_vector_dna<SampleDens,fixedIntWidth>::operator!=(const enc_vector_dna<SampleDens,fixedIntWidth>& v)const
{
    return !(*this == v);
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
void enc_vector_dna<SampleDens,fixedIntWidth>::swap(enc_vector_dna<SampleDens,fixedIntWidth>& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        m_z.swap(v.m_z);					// swap compressed bit streams
        m_sample_vals_and_pointer.swap(v.m_sample_vals_and_pointer);
        std::swap(m_elements, v.m_elements);// swap the number of elements
    }
}



template<uint32_t SampleDens, uint8_t fixedIntWidth>
template<class Container>
void enc_vector_dna<SampleDens,fixedIntWidth>::init(const Container& c)
{
    // clear BitVectors
    m_z.resize(0);
    m_elements = 0;
//	m_inc_start.resize(0);
    m_sample_vals_and_pointer.resize(0);
    if (c.empty())  // if c is empty there is nothing to do...
        return;
    typename Container::const_iterator	it		 	= c.begin(), end = c.end();
    typename Container::value_type 		v1			= *it, v2, max_value=0, max_sample_value=0, x;
    size_type samples=0;
    size_type z_size = 0, z_size_new=0, blocknr=0, newblocknr=0, writeblock=0;
    // invariant: blocknr != 8
    for (size_type i=0, no_sample=0; it != end; ++it,++i, --no_sample) {
        v2 = *it;
        if (!no_sample) { // add a sample
            no_sample = SampleDens;
            if (max_sample_value < v2) max_sample_value = v2;
            ++samples;
        } else {
            if (max_value < v2-v1) max_value = v2 - v1;
            if (v2 == v1) {
                throw std::logic_error("enc_vector_dna cannot decode adjacent equal values!");
            }
            z_size_new = z_size + coder::encoding_length(v2-v1);
            newblocknr = (z_size_new/64)%9;
            if (newblocknr == 8 or (newblocknr==0 and blocknr==7)) {  //if we reach or pass the sumblock
                z_size_new += 64; // add the bits for the sumblock
            }
            z_size = z_size_new;
            blocknr = (z_size/64)%9;
        }
        v1=v2;
    }
    if (blocknr==0 and (z_size%64)==0) {
        //no additional bits to add
    } else {
        assert(blocknr!=8);
        z_size += (8-blocknr)*64;
        blocknr = (z_size/64)%9;
        assert(blocknr==8);
        z_size += (64-(z_size%64)); // fill last bits
        blocknr = (z_size/64)%9;
        assert(blocknr==0 and (z_size%64)==0);
    }

//std::cerr<<"Calculate delta"<<std::endl;
    {
//		int_vector<> delta_c( c.size()-samples, 0, sizeof(typename Container::value_type)*8 ); // Vector for difference encoding of c
        if (max_sample_value > z_size+1)
            m_sample_vals_and_pointer.set_int_width(bit_magic::l1BP(max_sample_value) + 1);
        else
            m_sample_vals_and_pointer.set_int_width(bit_magic::l1BP(z_size+1) + 1);
        m_sample_vals_and_pointer.resize(2*samples+2); // add 2 for last entry
//		int_vector<0>::iterator d_it = delta_c.begin();
//		int_vector<0>::iterator sv_it = m_sample_vals_and_pointer.begin();
        typename enc_vector_dna_trait<fixedIntWidth>::int_vector_type::iterator sv_it = m_sample_vals_and_pointer.begin();
        z_size = 0;
        size_type no_sample=0;
        blocknr = 0; newblocknr = 0;
        for (it = c.begin(); it != end; ++it, --no_sample) {
            v2 = *it;
            if (!no_sample) { // add a sample
                no_sample = SampleDens;
                *sv_it = v2; ++sv_it;
                *sv_it = z_size; ++sv_it;
            } else {
                z_size_new = z_size + coder::encoding_length(v2-v1);
                newblocknr = (z_size_new/64)%9;
                if (newblocknr == 8 or (newblocknr==0 and blocknr==7)) {
                    z_size_new+=64;
                }
                z_size = z_size_new;
                blocknr = (z_size/64)%9;
            }
            v1=v2;
        }
        *sv_it = 0; ++sv_it;        // initialize
        *sv_it = z_size+1; ++sv_it; // last entry

        m_z.bit_resize(z_size);
        z_size = 0;
        uint64_t* z_data = coder::raw_data(m_z);
        uint8_t offset = 0;
        no_sample = 0;
        uint64_t sumblock[8]= {0};
        uint64_t* sbp=sumblock;
        bool overflow[8]= {false,false,false,false,false,false,false,false};
        bool* ovp = overflow;
        blocknr = 0; newblocknr = 0, writeblock = 0;
        for (it = c.begin(); it != end; ++it, --no_sample) {
            v2 = *it;
            if (!no_sample) { // add a sample
                no_sample = SampleDens;
            } else { // invariant: (z_data/64)%9!=8
                x = v2 - v1;
                coder::encode(x, z_data, offset);
                z_size_new = z_size + coder::encoding_length(x);
                newblocknr = (z_size_new/64)%9;
                writeblock = ((z_size_new-1)/64)%9;

                if (newblocknr == 8) { // if the next integer will be written in the sumblock
                    if ((z_size_new%64) > 0) {// if we have already written bits in the sumblock
                        *(z_data+1) = *z_data; // copy these bits in the next word
                    }
                    uint64_t sumblockword=0;
                    // calculate sumblock
                    if (writeblock == 7) {
                        write_sumblock(x, writeblock, sbp, ovp);
                        x=0;
                    }
                    for (int i=7; i>=0; --i) {
                        sumblockword<<=8;
                        sumblockword+=sumblock[i];
                        sumblock[i]=0;
                        overflow[i]=false;
                    }
                    *z_data = sumblockword;
                    ++z_data;
                    z_size_new += 64;//add 64 because we have added the sumblock
                } else if (newblocknr==0 and blocknr==7) {//
                    if ((z_size_new%64) > 0) {
                        *(z_data+1) = *z_data; // copy bit from block 0 to block 1
                    }
                    *z_data = *(z_data-1);	// copy the content of the sumblock to block 0
                    uint64_t sumblockword=0;
                    // calculate sumblock
                    for (int i=7; i>=0; --i) {
                        sumblockword<<=8;
                        sumblockword+=sumblock[i];
                        sumblock[i]=0;
                        overflow[i]=false;
                    }
                    *(z_data-1) = sumblockword;
                    ++z_data;
                    z_size_new += 64;// add 64 because we have added the sumblock
                }
                z_size = z_size_new;
                blocknr = (z_size/64)%9;
                writeblock = ((z_size-1)/64)%9;

                // sumblock[blocknr] contains the sum of fibonacci encoded words that ends in this block
                // or 0 if this sum is greater than 255.
                write_sumblock(x, writeblock, sbp, ovp);
//				if(!overflow[blocknr]){
//					if( x > 255 ){
//						overflow[blocknr] = true; sumblock[blocknr] = 0;
//					}else{
//						sumblock[blocknr] += x;
//						if(sumblock[blocknr]>255){ overflow[blocknr] = true; sumblock[blocknr] = 0; }
//					}
//				}
            }
            v1=v2;
        }
        if (blocknr==0 and (z_size%64)==0) {
            //no additional bits to add
        } else { // write last sumblock
            assert(blocknr!=8);
            uint64_t sumblockword=0;
            // calculate sumblock
            for (int i=7; i>=0; --i) {
                sumblockword<<=8;
                sumblockword+=sumblock[i];
                sumblock[i]=0;
                overflow[i]=false;
            }
            z_data += (8-blocknr);
            *z_data = sumblockword;
        }
    }
//	delta_c.resize(0);
//std::cerr<<"Calc rank"<<std::endl;
//std::cerr<<"Calc select"<<std::endl;
//std::cerr<<"Finished "<<std::endl;,
    m_elements = c.size();
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
enc_vector_dna<>::size_type enc_vector_dna<SampleDens,fixedIntWidth>::serialize(std::ostream& out) const
{
    size_type written_bytes = 0;
    out.write((char*) &m_elements, sizeof(m_elements));
    written_bytes += sizeof(m_elements);
    written_bytes += m_z.serialize(out);
    written_bytes += m_sample_vals_and_pointer.serialize(out);
    return written_bytes;
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
void enc_vector_dna<SampleDens,fixedIntWidth>::load(std::istream& in)
{
    in.read((char*) &m_elements, sizeof(m_elements));
    m_z.load(in);
    m_sample_vals_and_pointer.load(in);
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
const typename enc_vector_dna<SampleDens,fixedIntWidth>::const_iterator enc_vector_dna<SampleDens,fixedIntWidth>::begin()const
{
    return const_iterator(this, 0);
}

template<uint32_t SampleDens, uint8_t fixedIntWidth>
const typename enc_vector_dna<SampleDens,fixedIntWidth>::const_iterator enc_vector_dna<SampleDens,fixedIntWidth>::end()const
{
    return const_iterator(this, this->m_elements);
}


} // end namespace sdsl

#endif
