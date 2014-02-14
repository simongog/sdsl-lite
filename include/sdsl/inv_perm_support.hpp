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
/*! \file inv_perm_support.hpp
    \brief inv_perm_support.hpp contains a class which adds access to the
           inverse of a permutation.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_INV_PERM_SUPPORT
#define INCLUDED_SDSL_INV_PERM_SUPPORT

#include "int_vector.hpp"
#include "iterators.hpp"
#include "bit_vectors.hpp"
#include "rank_support.hpp"

namespace sdsl
{

//! Class inv_perm_support adds access to the inverse of a permutation.
/*!
 * \tparam t_s    Sampling parameter of the inverse permutation.
 * \tparam t_bv   Type of the bitvector used to indicate back-pointers.
 * \tparam t_rank Type of rank_support to rank the indicator bitvector.
 *
 * This support class adds access to the inverse of a permutation in at
 * most \(t_s\) steps. It takes about \(1/t_s \log n\) space, where \(n\)
 * is the size of the supported permutation.
 *
 * \par References
 *      [1] J. Munro, R. Raman, V. Raman, S. Rao: ,,Succinct representation
 *          of permutations'', Proceedings of ICALP 2003
 */
template<uint64_t t_s=32, class t_bv=bit_vector, class t_rank=typename bit_vector::rank_1_type>
class inv_perm_support
{
    public:
        typedef int_vector<>                                   iv_type;
        typedef iv_type::size_type                             size_type;
        typedef iv_type::value_type                            value_type;
        typedef iv_type::difference_type                       difference_type;
        typedef random_access_const_iterator<inv_perm_support> const_iterator;
        typedef t_bv                                           bit_vector_type;
        typedef t_rank                                         rank_type;
    private:
        const iv_type*  m_v = nullptr;   // pointer to supported permutation
        iv_type         m_back_pointer;  // back pointers
        bit_vector_type m_marked;        // back pointer marking
        rank_type       m_rank_marked;   // rank support for back pointer marking
    public:

        inv_perm_support() {};

        inv_perm_support(const inv_perm_support& p) : m_v(p.m_v),
            m_back_pointer(p.m_back_pointer), m_marked(p.m_marked),
            m_rank_marked(p.m_rank_marked) {
            m_rank_marked.set_vector(&m_marked);
        }

        inv_perm_support(inv_perm_support&& p) {
            *this = std::move(p);
        }

        //! Constructor
        inv_perm_support(const iv_type* v) : m_v(v) {
            bit_vector marked = bit_vector(m_v->size(), 0);
            bit_vector done   = bit_vector(m_v->size(), 0);

            size_type max_back_pointer = 0;
            for (size_type i = 0; i < m_v->size(); ++i) {
                if (!done[i]) {
                    done[i] = 1;
                    size_type back_pointer=i, j = i, j_new=0;
                    uint64_t  steps = 0, all_steps = 0;
                    while ((j_new=(*m_v)[j]) != i) {
                        j = j_new;
                        done[j] = 1;
                        ++steps; ++all_steps;
                        if (t_s == steps) {
                            max_back_pointer = std::max(max_back_pointer, back_pointer);
                            marked[j] = 1;
                            steps = 0;
                            back_pointer = j;
                        }
                    }
                    if (all_steps > t_s) {
                        marked[i] = 1;
                        max_back_pointer = std::max(max_back_pointer, back_pointer);
                    }
                }
            }

            m_marked = t_bv(std::move(marked));
            util::init_support(m_rank_marked, &m_marked);

            done   = bit_vector(m_v->size(), 0);
            size_type n_bp = m_rank_marked(m_v->size());
            m_back_pointer = int_vector<>(n_bp, 0, bits::hi(max_back_pointer)+1);

            for (size_type i = 0; i < m_v->size(); ++i) {
                if (!done[i]) {
                    done[i] = 1;
                    size_type back_pointer = i, j = i, j_new=0;
                    uint64_t  steps = 0, all_steps = 0;
                    while ((j_new=(*m_v)[j]) != i) {
                        j = j_new;
                        done[j] = 1;
                        ++steps; ++all_steps;
                        if (t_s == steps) {
                            m_back_pointer[m_rank_marked(j)] = back_pointer;
                            steps = 0;
                            back_pointer = j;
                        }
                    }
                    if (all_steps > t_s) {
                        m_back_pointer[m_rank_marked(i)] = back_pointer;
                    }
                }
            }
        }

        //! Access operator
        value_type operator[](size_type i) const {
            size_type j = i, j_new=0;
            while ((j_new=(*m_v)[j]) != i) {
                if (m_marked[j]) {
                    j = m_back_pointer[m_rank_marked(j)];
                    while ((j_new=(*m_v)[j]) != i) j = j_new;
                } else {
                    j = j_new;
                }
            }
            return j;
        }

        size_type size() const {
            return nullptr == m_v ? 0 : m_v->size();
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }

        void set_vector(const iv_type* v) { m_v = v; }

        //! Assignment operation
        inv_perm_support& operator=(const inv_perm_support& p) {
            if (this != &p) {
                m_v            = p.m_v;
                m_back_pointer = p.m_back_pointer;
                m_marked       = p.m_marked;
                m_rank_marked  = p.m_rank_marked;
                m_rank_marked.set_vector(&m_marked);
            }
            return *this;
        }

        //! Assignment move operation
        inv_perm_support& operator=(inv_perm_support&& p) {
            if (this != &p) {
                m_v            = std::move(p.m_v);
                m_back_pointer = std::move(p.m_back_pointer);
                m_marked       = std::move(p.m_marked);
                m_rank_marked  = std::move(p.m_rank_marked);
                m_rank_marked.set_vector(&m_marked);
            }
            return *this;
        }

        //! Swap operation
        void swap(inv_perm_support& p) {
            if (this != &p) {
                m_back_pointer.swap(p.m_back_pointer);
                m_marked.swap(p.m_marked);
                util::swap_support(m_rank_marked, p.m_rank_marked, &m_marked, &(p.m_marked));
            }
        }

        //! Serialize into stream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_back_pointer.serialize(out, child, "back_pointer");
            written_bytes += m_marked.serialize(out, child, "marked");
            written_bytes += m_rank_marked.serialize(out, child, "rank_marked");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load sampling from disk
        void load(std::istream& in) {
            m_back_pointer.load(in);
            m_marked.load(in);
            m_rank_marked.load(in, &m_marked);
        }
};

} // end namespace sdsl

#endif
