/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
/*! \file algorithms_for_compressed_suffix_trees.hpp
    \brief algorithms_for_compressed_suffix_trees.hpp contains algorithms for compressed suffix trees.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_ALGORITHMS_FOR_COMPRESSED_SUFFIX_TREES
#define INCLUDED_SDSL_ALGORITHMS_FOR_COMPRESSED_SUFFIX_TREES

#include "int_vector.hpp" // for bit_vector
#include "sorted_stack_support.hpp" // for construct_supercartesian_tree_bp
#include "sorted_multi_stack_support.hpp" // for first_p_index_construction
#include "util.hpp"
#include <stack> // for calculate_supercartesian_tree_bp


namespace sdsl
{

namespace algorithm
{

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param vec Random access container for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$ \Order{n \cdot \log n } \f$ bits.
 */
template<class RandomAccessContainer>
void construct_supercartesian_tree_bp(const RandomAccessContainer& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanaced parantheses to 2 n bits
    util::set_zero_bits(bp);
    std::stack<typename RandomAccessContainer::value_type> vec_stack;

    size_type k=0;
    for (size_type i=0; i < vec.size(); ++i) {
        typename RandomAccessContainer::value_type l = vec[i];
        if (minimum) {
            while (vec_stack.size() > 0 and l < vec_stack.top()) {
                vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
            }
        } else {
            while (vec_stack.size() > 0 and l > vec_stack.top()) {
                vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
            }
        }
        vec_stack.push(l);
        bp[k++] = 1; // writing an opening  parenthesis
    }
    while (vec_stack.size() > 0) {
        vec_stack.pop();
        bp[k++] = 0; // writing a closing parenthesis
    }
    assert(k == 2*vec.size());
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param vec Random access container for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{n}\f$ bits, by the stack_support described in the paper "Optimal Succinctness For Range Minimum Queries" of Johannes Fischer.
 */
// TODO: sorted_multi_stack_support einbauen, RandomAccessContainer durch int_vector_file_buffer ersetzen
template<class RandomAccessContainer>
void construct_supercartesian_tree_bp_succinct(const RandomAccessContainer& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanced parentheses to 2 n bits
	if ( vec.size() > 0 ) {
		util::set_zero_bits(bp);
		sorted_stack_support vec_stack(vec.size()); // <- das ist ein Problem fuer int_vector_file_buffer

		size_type k=0;
		if ( minimum ) {
			bp[k++] = 1;
			for (size_type i=1; i < vec.size(); ++i) {
				if (vec[i] < vec[i-1]) {
					++k;
					while (vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()]) {
						vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zero
					}
				} else {
					vec_stack.push(i-1); // "lazy stack" trick: speed-up ca. 25%
				}
				bp[k++] = 1; // writing an opening  parenthesis
			}
			/*
			vec_stack.push(0);
			bp[k++] = 1;
			for(size_type i=1,j, start_run=1; i < vec.size(); ++i){
				if( vec[i] < vec[i-1] ){
					j = i;
					while( --j >= start_run and vec[i] < vec[j]) ++k;
					while(start_run <= j){	// auf den stack pushen
						vec_stack.push(start_run++);
					}
					while( vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()] ){
						vec_stack.pop(); ++k;
					}
					start_run = i;
				}
				bp[k++] = 1;
			}
			*/
		} else {
			// hier noch ohne "lazy stack" trick
			for (size_type i=0; i < vec.size(); ++i) {
				while (vec_stack.size() > 0 and vec[i] > vec[vec_stack.top()]) {
					vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
				}
				vec_stack.push(i);
				bp[k++] = 1; // writing an opening  parenthesis
			}
		}
#ifdef SDSL_DEBUG
	// not necessary as bp is already initialized to zero
		while (!vec_stack.empty()) {
			vec_stack.pop();
			bp[k++] = 0; // writing a closing parenthesis
		}
		assert(k == 2*vec.size());
#endif
	}
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param lcp_buf int_vector_file_buffer of the LCP Array for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{2n}\f$ bits, by the multi_stack_support
 */
template<uint8_t fixedIntWidth>
void construct_supercartesian_tree_bp_succinct(int_vector_file_buffer<fixedIntWidth>& lcp_buf, bit_vector& bp, const bool minimum=true)
{
    typedef int_vector_size_type size_type;
    lcp_buf.reset();
    size_type n = lcp_buf.int_vector_size;
    bp.resize(2*n);      // resize bit vector for balanced parentheses to 2 n bits
    if (n == 0)	// if n == 0 we are done
        return;
    util::set_zero_bits(bp);
    sorted_multi_stack_support vec_stack(n);

    size_type k=0;
    if (minimum) {
        bp[k++] = 1;
        size_type r = lcp_buf.load_next_block();
        size_type last = lcp_buf[0];
        for (size_type i=1, r_sum = 0, x; r_sum < n;) {
            for (; i < r_sum +r; ++i) {
                x = lcp_buf[i-r_sum];
                if (x < last) {
                    ++k; // writing a closing parenthesis for last
                    while (!vec_stack.empty() and x < vec_stack.top()) {
                        vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zeros
                    }
                } else {
                    vec_stack.push(last); // "lazy stack" trick: Beschleunigung: ca 25 %
                }
                bp[k++] = 1; // writing an opening parenthesis
                last = x;
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
    } else {
        // hier noch ohne "lazy stack" trick
        for (size_type i=0, r_sum = 0, r = lcp_buf.load_next_block(), x; r_sum < n;) {
            for (; i < r_sum +r; ++i) {
                x = lcp_buf[i-r_sum];
                while (!vec_stack.empty() and x > vec_stack.top()) {
                    vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zeros
                }
                vec_stack.push(x);
                bp[k++] = 1; // writing an opening parenthesis
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
    }
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009) and the first_child bit_vector
/*! \param lcp_buf int_vector_file_buffer for the lcp array for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param bp_fc Reference to the first child bit_vector of bp.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{2n}\f$ bits, by the multi_stack_support
 */
template<uint8_t fixedIntWidth>
int_vector_size_type construct_supercartesian_tree_bp_succinct_and_first_child(int_vector_file_buffer<fixedIntWidth>& lcp_buf, bit_vector& bp, bit_vector& bp_fc, const bool minimum=true)
{
    typedef int_vector_size_type size_type;
    lcp_buf.reset();
    size_type n = lcp_buf.int_vector_size;
    bp.resize(2*n);      // resize bit vector for balanaced parantheses to 2 n bits
    bp_fc.resize(n);
    if (n == 0)	// if n == 0 we are done
        return 0;
    size_type fc_cnt=0; // first child counter
    util::set_zero_bits(bp);
    util::set_zero_bits(bp_fc);
    sorted_multi_stack_support vec_stack(n);

    size_type k=0;
    size_type k_fc=0; // first child index
    if (minimum) {
        // hier noch ohne "lazy stack" trick
        for (size_type i=0, r_sum = 0, r = lcp_buf.load_next_block(), x; r_sum < n;) {
            for (; i < r_sum +r; ++i) {
                x = lcp_buf[i-r_sum];
                while (!vec_stack.empty() and x < vec_stack.top()) {
                    if (vec_stack.pop()) {
                        bp_fc[k_fc] = 1;
                        ++fc_cnt;
                    }
                    ++k; // writing a closing parenthesis, bp is already initialized to zeros
                    ++k_fc; // write a bit in first_child
                }
                vec_stack.push(x);
                bp[k++] = 1; // writing an opening parenthesis
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }

    } else {
        // hier noch ohne "lazy stack" trick
        for (size_type i=0, r_sum = 0, r = lcp_buf.load_next_block(), x; r_sum < n;) {
            for (; i < r_sum +r; ++i) {
                x = lcp_buf[i-r_sum];
                while (!vec_stack.empty() and x > vec_stack.top()) {
                    if (vec_stack.pop()) {
                        bp_fc[k_fc] = 1;
                        ++fc_cnt;
                    }
                    ++k; // writing a closing parenthesis, bp is already initialized to zeros
                    ++k_fc; // write a bit in first_child
                }
                vec_stack.push(x);
                bp[k++] = 1; // writing an opening parenthesis
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
    }
    while (!vec_stack.empty()) {
        if (vec_stack.pop()) {
            bp_fc[k_fc] = 1;
            ++fc_cnt;
        }
        // writing a closing parenthesis in bp, not necessary as bp is initalized with zeros
        ++k;
        ++k_fc;
    }
//	assert( k == 2*vec.size() );
    return fc_cnt;
}

template<uint8_t fixedIntWidth>
int_vector_size_type construct_first_child_lcp(int_vector_file_buffer<fixedIntWidth>& lcp_buf, int_vector<fixedIntWidth>& fc_lcp, int_vector_size_type first_child_size=0)
{
    typedef int_vector_size_type size_type;
    lcp_buf.reset();
    size_type n = lcp_buf.int_vector_size;
    if (n == 0) {	// if n == 0 we are done
        fc_lcp = int_vector<>(0);
        return 0;
    }
    if (first_child_size == 0) {
        util::assign(fc_lcp, int_vector<>(n, 0, bit_magic::l1BP(n)+1));
    }

    size_type fc_cnt=0; // first child counter
    sorted_multi_stack_support vec_stack(n);
    size_type y;
    for (size_type i=0, r_sum = 0, r = lcp_buf.load_next_block(), x; r_sum < n;) {
        for (; i < r_sum +r; ++i) {
            x = lcp_buf[i-r_sum];
            while (!vec_stack.empty() and x < vec_stack.top()) {
                y = vec_stack.top();
                if (vec_stack.pop()) {
                    fc_lcp[fc_cnt++] = y;
                }
            }
            vec_stack.push(x);
        }
        r_sum += r; r = lcp_buf.load_next_block();
    }

    while (!vec_stack.empty()) {
        y = vec_stack.top();
        if (vec_stack.pop()) {
            fc_lcp[fc_cnt++] = y;
        }
        // writing a closing parenthesis in bp, not necessary as bp is initalized with zeros
    }
    if (fc_cnt < fc_lcp.size()) {
        fc_lcp.resize(fc_cnt);
    }
    return fc_cnt;
}


template<uint32_t SampleDens, uint8_t fixedIntWidth, uint8_t bwtIntWidth>
int_vector_size_type construct_first_child_and_lf_lcp(int_vector_file_buffer<fixedIntWidth>& lcp_buf,
        int_vector_file_buffer<bwtIntWidth>& bwt_buf,
        std::string small_lcp_file_name,
        std::string big_lcp_file_name,
        int_vector<>& big_lcp,
        int_vector_size_type first_child_size=0)
{
    typedef int_vector_size_type size_type;
    const size_type M = 255;	// limit for values represented in the small LCP part
    size_type 		buf_len = 1000000;
    lcp_buf.reset(buf_len);
    bwt_buf.reset(buf_len);
    size_type n = lcp_buf.int_vector_size;

    std::ofstream sml_lcp_out(small_lcp_file_name.c_str());
    size_type bit_size = 8*n;
    sml_lcp_out.write((char*) &bit_size, sizeof(bit_size));

    std::ofstream big_lcp_out(big_lcp_file_name.c_str());

    size_type fc_cnt = 0; // number of lcp values at the first child r
    size_type fc_cnt_big = 0; // number of lcp values at the first child which are big and not reducable
    size_type fc_cnt_big2 = 0;
    sorted_multi_stack_support vec_stack(n); // occupies 2n bits
    bit_vector is_big_and_not_reducable(n); // initialized with 0s
    bool is_one_big_and_not_reducable = false; // all positions have to be reducable

    size_type y, max_lcp=0;
    uint64_t last_bwti=0, val;
    for (size_type i=0, r_sum = 0, r = 0, x; r_sum < n;) {
        for (; i < r_sum + r; ++i) {
            x = lcp_buf[i-r_sum];
            is_one_big_and_not_reducable = false;

            while (!vec_stack.empty() and x < vec_stack.top()) {
                y = vec_stack.top();
                is_one_big_and_not_reducable |= is_big_and_not_reducable[vec_stack.size()-1];
                if (vec_stack.pop()) { // if y was the last copy of y on the stack
                    if (y > M-2) {
                        if (is_one_big_and_not_reducable) {
                            val = M;
                            big_lcp_out.write((char*)&y, sizeof(y));
                            ++fc_cnt_big;
                            if (y > max_lcp) max_lcp = y;
                        } else {
                            val = M-1;
                            ++fc_cnt_big2;
                        }
                    } else {
                        val = y;
                    }
                    sml_lcp_out.write((const char*)&val, 1);
                    ++fc_cnt;
                    is_one_big_and_not_reducable = false;
                }
            }
            if (x > M-2 and (0 == i or last_bwti != bwt_buf[i - r_sum] or x % SampleDens == 0)) {
                is_big_and_not_reducable[vec_stack.size()] = 1;
            } else {
                is_big_and_not_reducable[vec_stack.size()] = 0;
            }
            vec_stack.push(x);
            last_bwti = bwt_buf[i - r_sum];
        }
        r_sum += r; r = lcp_buf.load_next_block();
        bwt_buf.load_next_block();
    }

    while (!vec_stack.empty()) {
        y = vec_stack.top();
        if (vec_stack.pop()) {
            if (y > M-2) {
                if (is_big_and_not_reducable[vec_stack.size()]) {
                    val = M;
                    big_lcp_out.write((char*)&y, sizeof(y));
                    ++fc_cnt_big;
                    if (y > max_lcp) max_lcp = y;
                } else {
                    val = M-1;
                    ++fc_cnt_big2;
                }
            } else {
                val = y;
            }
            sml_lcp_out.write((const char*)&val, 1);
            ++fc_cnt;
        }
    }

    // write number of elements of sml_lcp into the out file stream
    sml_lcp_out.seekp(0);
    bit_size = 8*fc_cnt;
    sml_lcp_out.write((char*) &bit_size, sizeof(bit_size));
    sml_lcp_out.close();

    big_lcp_out.close();
    std::ifstream big_lcp_in(big_lcp_file_name.c_str());
    big_lcp.set_int_width(bit_magic::l1BP(max_lcp)+1);
    big_lcp.resize(fc_cnt_big);

    for (size_type i=0; i<fc_cnt_big; ++i) {
        big_lcp_in.read((char*)&y, sizeof(y));
        big_lcp[i] = y;
    }
    big_lcp_in.close();
    std::remove(big_lcp_file_name.c_str());

//    std::cout<<"#n="<<n<<" fc_cnt="<<fc_cnt<<" fc_cnt_big="<<fc_cnt_big<<" fc_cnt_big2="<<fc_cnt_big2<<std::endl;

    return fc_cnt;
}



template<class RandomAccessContainer>
void construct_supercartesian_tree_bp_succinct2(const RandomAccessContainer& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanaced parantheses to 2 n bits
    util::set_zero_bits(bp);
    sorted_stack_support vec_stack(vec.size()); // <- das ist ein Problem fuer int_vector_file_buffer

    size_type k=0;
//	uint64_t wbuf=0;
    for (size_type i=0/*, cnt64=0*/; i < vec.size(); ++i) {
        while (vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()]) {
            vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
        }
        vec_stack.push(i);
        bp[k++] = 1; // writing an opening  parenthesis
        while (i+1 < vec.size() and vec[i+1] >= vec[i]) {
            vec_stack.push(++i);
            bp[k++];
        }
    }
#ifdef SDSL_DEBUG
// not neccessary as bp is already initialized to zero
    while (vec_stack.size() > 0) {
        vec_stack.pop();
        bp[k++] = 0; // writing a closing parenthesis
    }
    assert(k == 2*vec.size());
#endif
}

template<class RandomAccessContainer>
typename RandomAccessContainer::size_type construct_first_p_index(const RandomAccessContainer& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    size_type nr_of_first_indices = 0;
    bp = bit_vector(vec.size(), 0);
//	std::cerr<<"bp.size()="<<bp.size()<<std::endl;
    sorted_stack_support vec_stack(vec.size());
    size_type k=0;
    for (size_type i=0, t; i < vec.size(); ++i) {
        if (minimum) {
            while (vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()]) {
                t = vec[vec_stack.top()];
                vec_stack.pop();
                if (vec_stack.size() == 0 or t != vec[vec_stack.top()]) {
                    bp[k] = 1;
                    ++nr_of_first_indices;
                }
                ++k;

            }
        } else {
            while (vec_stack.size() > 0 and vec[i] > vec[vec_stack.top()]) {
                t = vec[vec_stack.top()];
                vec_stack.pop();
                if (vec_stack.size() == 0 or t != vec[vec_stack.top()]) {
                    bp[k] = 1;
                    ++nr_of_first_indices;
                }
                ++k;
            }
        }
        vec_stack.push(i);
    }
    while (vec_stack.size() > 0) {
        size_type t = vec[vec_stack.top()];
        vec_stack.pop();
        if (vec_stack.size() == 0 or t != vec[vec_stack.top()]) {
            bp[k] = 1;
            ++nr_of_first_indices;
        }
        ++k;
    }
    assert(k == vec.size());
    return nr_of_first_indices;
}

template<uint8_t fixedIntWidth>
bit_vector::size_type construct_first_p_index(int_vector_file_buffer<fixedIntWidth>& lcp_buf, bit_vector& bp, const bool minimum=true)
{
    typedef bit_vector::size_type size_type;
    size_type nr_of_first_indices = 0;
    lcp_buf.reset();
    size_type n = lcp_buf.int_vector_size;

    bp = bit_vector(n, 0);
    sorted_multi_stack_support vec_stack(n);
    size_type k=0;

    if (minimum) {
        for (size_type i = 0, r_sum = 0, r = lcp_buf.load_next_block(),x; r_sum < n;) {
            for (; i<r_sum+r; ++i) {
                x = lcp_buf[i-r_sum];
                while (!vec_stack.empty() and x < vec_stack.top()) {
                    if (vec_stack.pop()) {
                        bp[k] = 1;
                        ++nr_of_first_indices;
                    }
                    ++k;
                }
                vec_stack.push(x);
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
    } else {
        for (size_type i = 0, r_sum = 0, r = lcp_buf.load_next_block(),x; r_sum < n;) {
            for (; i<r_sum+r; ++i) {
                x = lcp_buf[i-r_sum];
                while (!vec_stack.empty() and x > vec_stack.top()) {
                    if (vec_stack.pop()) {
                        bp[k] = 1;
                        ++nr_of_first_indices;
                    }
                    ++k;
                }
                vec_stack.push(x);
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
    }

    while (!vec_stack.empty()) {
        if (vec_stack.pop()) {
            bp[k] = 1;
            ++nr_of_first_indices;
        }
        ++k;
    }
//	assert( k == vec.size() );
    return nr_of_first_indices;
}

}// end namespace algorithm

}// end namespace sdsl

#endif

