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
/*! \file wm_int.hpp
    \brief wm_int.hpp contains a specialized class for a wavelet tree for
           sequences over large alphabets.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_WM_INT
#define INCLUDED_SDSL_WM_INT

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "wt_helper.hpp"
#include "util.hpp"
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A wavelet tree class for integer sequences.
/*!
 * \tparam t_bitvector   Type of the bitvector used for representing the wavelet tree.
 * \tparam t_rank        Type of the support structure for rank on pattern `1`.
 * \tparam t_select      Type of the support structure for select on pattern `1`.
 * \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 * This wavelet tree variant does not store the two children of a node v aligned
 * with v; it is also known as wavelet matrix.
 *
 * \par References
 *      [1] F. Claude, G. Navarro: ,,The Wavelet Matrix'', Proceedings of
 *          SPIRE 2012.
 *
 *   @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
class wm_int
{
    public:

        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef typename t_bitvector::difference_type difference_type;
        typedef random_access_const_iterator<wm_int> const_iterator;
        typedef const_iterator                       iterator;
        typedef t_bitvector                          bit_vector_type;
        typedef t_rank                               rank_1_type;
        typedef t_select                             select_1_type;
        typedef t_select_zero                        select_0_type;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        enum 	{lex_ordered=0};

        typedef std::pair<value_type, size_type>     point_type;
        typedef std::vector<point_type>              point_vec_type;
        typedef std::pair<size_type, point_vec_type> r2d_res_type;

        struct node_type;


    protected:

        size_type              m_size  = 0;
        size_type              m_sigma = 0;    //<- \f$ |\Sigma| \f$
        bit_vector_type        m_tree;         // bit vector to store the wavelet tree
        rank_1_type            m_tree_rank;    // rank support for the wavelet tree bit vector
        select_1_type          m_tree_select1; // select support for the wavelet tree bit vector
        select_0_type          m_tree_select0;
        uint32_t               m_max_level = 0;
        int_vector<64>         m_zero_cnt;     // m_zero_cnt[i] contains the number of zeros in level i
        int_vector<64>         m_rank_level;   // m_rank_level[i] contains m_tree_rank(i*size())

        void copy(const wm_int& wt)
        {
            m_size          = wt.m_size;
            m_sigma         = wt.m_sigma;
            m_tree          = wt.m_tree;
            m_tree_rank     = wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1  = wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0  = wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            m_max_level     = wt.m_max_level;
            m_zero_cnt      = wt.m_zero_cnt;
            m_rank_level    = wt.m_rank_level;
        }

    public:

        const size_type&       sigma = m_sigma;         //!< Effective alphabet size of the wavelet tree.
        const bit_vector_type& tree  = m_tree;          //!< A concatenation of all bit vectors of the wavelet tree.
        const uint32_t&        max_level = m_max_level; //!< Maximal level of the wavelet tree.

        //! Default constructor
        wm_int()
        {
        };

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wm_int should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         *  \param max_level   Maximal level of the wavelet tree. If set to 0, determined automatically.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *    \par Space complexity
         *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wm_int(int_vector_buffer<int_width>& buf, size_type size,
               uint32_t max_level=0) : m_size(size)
        {
            if (0 == m_size)
                return;
            size_type n = buf.size();  // set n
            if (n < m_size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
                return;
            }
            m_sigma = 0; // init sigma

            int_vector<int_width> rac(m_size, 0, buf.width());  // initialize rac

            value_type x = 1;  // variable for the biggest value in rac
            for (size_type i=0; i < m_size; ++i) { // detect the largest value in rac
                if (buf[i] > x)
                    x = buf[i];
                rac[i] = buf[i];
            }

            if (max_level == 0) {
                m_max_level = bits::hi(x)+1; // we need max_level bits to represent all values in the range [0..x]
            } else {
                m_max_level = max_level;
            }


            std::string tree_out_buf_file_name = tmp_file(buf.filename(), "_m_tree");
            osfstream tree_out_buf(tree_out_buf_file_name, std::ios::binary | std::ios::trunc | std::ios::out);   // open buffer for tree
            size_type bit_size = m_size*m_max_level;
            tree_out_buf.write((char*) &bit_size, sizeof(bit_size));    // write size of bit_vector

            std::string zero_buf_file_name = tmp_file(buf.filename(), "_zero_buf");

            size_type tree_pos = 0;
            uint64_t tree_word = 0;

            m_zero_cnt = int_vector<64>(m_max_level, 0); // zeros at level i

            for (uint32_t k=0; k<m_max_level; ++k) {
                uint8_t        width = m_max_level-k-1;
                const uint64_t mask  = 1ULL<<width;
                uint64_t       x     = 0;
                size_type      zeros = 0;
                int_vector_buffer<> zero_buf(zero_buf_file_name, std::ios::out, 1024*1024, m_max_level);
                for (size_t i=0; i<m_size; ++i) {
                    x = rac[i];
                    if (x&mask) {
                        tree_word |= (1ULL << (tree_pos&0x3FULL));
                        zero_buf.push_back(x);
                    } else {
                        rac[zeros++ ] = x;
                    }
                    ++tree_pos;
                    if ((tree_pos & 0x3FULL) == 0) { // if tree_pos % 64 == 0 write old word
                        tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
                        tree_word = 0;
                    }
                }
                m_zero_cnt[k] = zeros;
                for (size_t i=zeros; i<m_size; ++i) {
                    rac[i] = zero_buf[i-zeros];
                }
            }
            if ((tree_pos & 0x3FULL) != 0) { // if tree_pos % 64 > 0 => there are remaining entries we have to write
                tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
            }
            sdsl::remove(zero_buf_file_name);
            tree_out_buf.close();
            m_sigma = std::unique(rac.begin(), rac.end()) - rac.begin();
            rac.resize(0);
            bit_vector tree;
            load_from_file(tree, tree_out_buf_file_name);
            sdsl::remove(tree_out_buf_file_name);
            m_tree = bit_vector_type(std::move(tree));
            util::init_support(m_tree_rank, &m_tree);
            util::init_support(m_tree_select0, &m_tree);
            util::init_support(m_tree_select1, &m_tree);
            m_rank_level = int_vector<64>(m_max_level, 0);
            for (uint32_t k=0; k<m_rank_level.size(); ++k) {
                m_rank_level[k] = m_tree_rank(k*m_size);
            }
        }

        //! Copy constructor
        wm_int(const wm_int& wt)
        {
            copy(wt);
        }

        //! Copy constructor
        wm_int(wm_int&& wt)
        {
            *this = std::move(wt);
        }

        //! Assignment operator
        wm_int& operator=(const wm_int& wt)
        {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        wm_int& operator=(wm_int&& wt)
        {
            if (this != &wt) {
                m_size          = wt.m_size;
                m_sigma         = wt.m_sigma;
                m_tree          = std::move(wt.m_tree);
                m_tree_rank     = std::move(wt.m_tree_rank);
                m_tree_rank.set_vector(&m_tree);
                m_tree_select1  = std::move(wt.m_tree_select1);
                m_tree_select1.set_vector(&m_tree);
                m_tree_select0  = std::move(wt.m_tree_select0);
                m_tree_select0.set_vector(&m_tree);
                m_max_level     = std::move(wt.m_max_level);
                m_zero_cnt      = std::move(wt.m_zero_cnt);
                m_rank_level    = std::move(wt.m_rank_level);
            }
            return *this;
        }

        //! Swap operator
        void swap(wm_int& wt)
        {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
                util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree));
                std::swap(m_max_level,  wt.m_max_level);
                m_zero_cnt.swap(wt.m_zero_cnt);
                m_rank_level.swap(wt.m_rank_level);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const
        {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const
        {
            return m_size == 0;
        }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const
        {
            assert(i < size());
            value_type res = 0;
            for (uint32_t k=0; k < m_max_level; ++k) {
                res <<= 1;
                size_type rank_ones = m_tree_rank(i) - m_rank_level[k];
                if (m_tree[i]) { // one at position i => follow right child
                    i = (k+1)*m_size + m_zero_cnt[k] + rank_ones;
                    res |= 1;
                } else { // zero at position i => follow left child
                    auto rank_zeros = (i - k*m_size) - rank_ones;
                    i = (k+1)*m_size + rank_zeros;
                }
            }
            return res;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const
        {
            assert(i <= size());
            if (((1ULL)<<(m_max_level))<=c) { // c is greater than any symbol in wt
                return 0;
            }
            size_type b = 0; // start position of the interval
            uint64_t mask = (1ULL) << (m_max_level-1);
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + i) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [level_b..b)
                if (c & mask) { // search for a one at this level
                    i = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                } else { // search for a zero at this level
                    i = i-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
            }
            return i;
        };

        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(wt[i],i),wt[i])
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type>
        inverse_select(size_type i)const
        {
            assert(i < size());
            value_type c = 0;
            size_type b = 0; // start position of the interval
            uint64_t mask = (1ULL) << (m_max_level-1);
            for (uint32_t k=0; k < m_max_level; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + i) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [level_b..b)
                c<<=1;
                if (m_tree[b+i]) { // go to the right child
                    i = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                    c|=1;
                } else { // go to the left child
                    i = i-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
            }
            return std::make_pair(i,c);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const
        {
            assert(1 <= i and i <= rank(size(), c));
            uint64_t mask = 1ULL << (m_max_level-1);
            int_vector<64> m_path_off(max_level+1);
            int_vector<64> m_path_rank_off(max_level+1);
            m_path_off[0] = m_path_rank_off[0] = 0;
            size_type b = 0; // start position of the interval
            size_type r = i;
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + r) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [0..b)
                if (c & mask) { // search for a one at this level
                    r = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                } else { // search for a zero at this level
                    r = r-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
                m_path_off[k+1] = b;
                m_path_rank_off[k] = rank_b;
            }
            mask = 1ULL;
            for (uint32_t k=m_max_level; k>0; --k) {
                b = m_path_off[k-1];
                size_type rank_b = m_path_rank_off[k-1];
                if (c & mask) { // right child => search i'th one
                    i = m_tree_select1(rank_b + i) - b + 1;
                } else { // left child => search i'th zero
                    i = m_tree_select0(b - rank_b + i) - b + 1;
                }
                mask <<= 1;
            }
            return i-1;
        };

        //! range_search_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
        /*! \param lb     Left bound of index interval (inclusive)
         *  \param rb     Right bound of index interval (inclusive)
         *  \param vlb    Left bound of value interval (inclusive)
         *  \param vrb    Right bound of value interval (inclusive)
         *  \param report Should the matching points be returned?
         *  \return Pair (#of found points, vector of points), the vector is empty when
         *          report = false.
         */
        std::pair<size_type, std::vector<std::pair<value_type, size_type>>>
        range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                        bool report=true) const
        {

            if (vrb > (1ULL << m_max_level))
                vrb = (1ULL << m_max_level);
            if (vlb > vrb)
                return make_pair(0, point_vec_type());
            size_type cnt_answers = 0;
            point_vec_type point_vec;
            if (lb <= rb) {
                std::vector<size_type> is(m_max_level+1);
                std::vector<size_type> rank_off(m_max_level+1);
                _range_search_2d(root(), {lb, rb}, vlb, vrb, 0, is,
                                 rank_off, point_vec, report, cnt_answers);
            }
            return make_pair(cnt_answers, point_vec);
        }

        void
        _range_search_2d(node_type v, range_type r, value_type vlb,
                         value_type vrb, size_type ilb, std::vector<size_type>& is,
                         std::vector<size_type>& rank_off, point_vec_type& point_vec,
                         bool report, size_type& cnt_answers)
        const
        {
            using std::get;
            if (get<0>(r) > get<1>(r))
                return;
            is[v.level] = v.offset + get<0>(r);

            if (v.level == m_max_level) {
                for (size_type j=1; j <= sdsl::size(r) and report; ++j) {
                    size_type i = j;
                    size_type c = v.sym;
                    for (uint32_t k=m_max_level; k>0; --k) {
                        size_type offset = is[k-1];
                        size_type rank_offset = rank_off[k-1];
                        if (c&1) {
                            i = m_tree_select1(rank_offset+i)-offset+1;
                        } else {
                            i = m_tree_select0(offset-rank_offset+i)-offset+1;
                        }
                        c >>= 1;
                    }
                    point_vec.emplace_back(is[0]+i-1, v.sym);
                }
                cnt_answers += sdsl::size(r);
                return;
            } else {
                rank_off[v.level] = m_tree_rank(is[v.level]);
            }
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and  vlb < mid and mid) {
                _range_search_2d(get<0>(c_v),get<0>(c_r), vlb,
                                 std::min(vrb,mid-1), ilb, is, rank_off,
                                 point_vec, report, cnt_answers);
            }
            if (!sdsl::empty(get<1>(c_r)) and vrb >= mid) {
                _range_search_2d(get<1>(c_v), get<1>(c_r), std::max(mid, vlb),
                                 vrb, mid, is, rank_off, point_vec, report,
                                 cnt_answers);
            }
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const
        {
            return const_iterator(this, size());
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += m_tree.serialize(out, child, "tree");
            written_bytes += m_tree_rank.serialize(out, child, "tree_rank");
            written_bytes += m_tree_select1.serialize(out, child, "tree_select_1");
            written_bytes += m_tree_select0.serialize(out, child, "tree_select_0");
            written_bytes += write_member(m_max_level, out, child, "max_level");
            written_bytes += m_zero_cnt.serialize(out, child, "zero_cnt");
            written_bytes += m_rank_level.serialize(out, child, "rank_level");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            read_member(m_sigma, in);
            m_tree.load(in);
            m_tree_rank.load(in, &m_tree);
            m_tree_select1.load(in, &m_tree);
            m_tree_select0.load(in, &m_tree);
            read_member(m_max_level, in);
            m_zero_cnt.load(in);
            m_rank_level.load(in);
        }

        //! Represents a node in the wavelet tree
        struct node_type {
            size_type  offset   = 0;
            size_type  size     = 0;
            size_type  level    = 0;
            value_type sym      = 0;

            // Default constructor
            node_type(size_type o=0, size_type sz=0, size_type l=0,
                      value_type sy=0) :
                offset(o), size(sz), level(l), sym(sy) {}

            // Copy constructor
            node_type(const node_type&) = default;

            // Move copy constructor
            node_type(node_type&&) = default;

            // Assignment operator
            node_type& operator=(const node_type&) = default;

            // Move assignment operator
            node_type& operator=(node_type&&) = default;

            // Comparator operator
            bool operator==(const node_type& v) const
            {
                return offset == v.offset;
            }

            // Smaller operator
            bool operator<(const node_type& v) const
            {
                return offset < v.offset;
            }

            // Greater operator
            bool operator>(const node_type& v) const
            {
                return offset > v.offset;
            }
        };


        //! Checks if the node is a leaf node
        bool is_leaf(const node_type& v) const
        {
            return v.level == m_max_level;
        }

        //! Symbol of leaf node v
        value_type sym(const node_type& v) const
        {
            return v.sym;
        }

        //! Random access container to bitvector of node v
        auto bit_vec(const node_type& v) const -> node_bv_container<t_bitvector> {
            return node_bv_container<t_bitvector>(begin(v), end(v));
        }

        //! Random access container to sequence of node v
        auto seq(const node_type& v) const -> random_access_container<std::function<value_type(size_type)>> {
            return random_access_container<std::function<value_type(size_type)>>([&v, this](size_type i)
            {
                node_type vv = v;
                while (!is_leaf(vv)) {
                    auto vs = expand(vv);
                    auto rs = expand(vv, {0, i});
                    bool bit = *(begin(vv)+i);
                    i = std::get<1>(rs[bit]);
                    vv = vs[bit];
                }
                return sym(vv);
            }, size(v));
        }

        //! Indicates if node v is empty
        bool empty(const node_type& v) const
        {
            return v.size == (size_type)0;
        }

        //! Return the size of node v
        auto size(const node_type& v) const -> decltype(v.size)
        {
            return v.size;
        }

        //! Return the root node
        node_type root() const
        {
            return node_type(0, m_size, 0, 0);
        }

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::array<node_type, 2>
        expand(const node_type& v) const
        {
            node_type v_right = v;
            return expand(std::move(v_right));
        }

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::array<node_type, 2>
        expand(node_type&& v) const
        {
            node_type v_left;
            size_type rank_b = m_tree_rank(v.offset);
            size_type ones   = m_tree_rank(v.offset+v.size)-rank_b; // ones in [b..size)
            size_type ones_p = rank_b - m_rank_level[v.level];      // ones in [level_b..b)

            v_left.offset = (v.level+1)*m_size + (v.offset - v.level*m_size) - ones_p;
            v_left.size   = v.size - ones;
            v_left.level  = v.level + 1;
            v_left.sym    = v.sym<<1;

            v.offset = (v.level+1)*m_size + m_zero_cnt[v.level] + ones_p;
            v.size   = ones;
            v.level  = v.level + 1;
            v.sym    = (v.sym<<1)|1;

            return {std::move(v_left), v};
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_vec_type, 2>
        expand(const node_type& v,
               const range_vec_type& ranges) const
        {
            auto ranges_copy = ranges;
            return expand(v, std::move(ranges_copy));
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_vec_type, 2>
        expand(const node_type& v,
               range_vec_type&& ranges) const
        {
            auto v_sp_rank = m_tree_rank(v.offset);  // this is already calculated in expand(v)
            range_vec_type res(ranges.size());
            size_t i = 0;
            for (auto& r : ranges) {
                auto sp_rank    = m_tree_rank(v.offset + r[0]);
                auto right_size = m_tree_rank(v.offset + r[1] + 1)
                                  - sp_rank;
                auto left_size  = (r[1]-r[0]+1)-right_size;

                auto right_sp = sp_rank - v_sp_rank;
                auto left_sp  = r[0] - right_sp;

                r = {left_sp, left_sp + left_size - 1};
                res[i++] = {right_sp, right_sp + right_size - 1};
            }
            return {ranges, std::move(res)};
        }

        //! Returns for a range its left and right child ranges
        /*! \param v An inner node of an wavelet tree.
         *  \param r A ranges [s,e], such that [s,e] is
         *           contained in v=[v_s,v_e].
         *  \return A range pair. The first element of the
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_type, 2>
        expand(const node_type& v, const range_type& r) const
        {
            auto v_sp_rank = m_tree_rank(v.offset);  // this is already calculated in expand(v)
            auto sp_rank    = m_tree_rank(v.offset + r[0]);
            auto right_size = m_tree_rank(v.offset + r[1] + 1)
                              - sp_rank;
            auto left_size  = (r[1]-r[0]+1)-right_size;

            auto right_sp = sp_rank - v_sp_rank;
            auto left_sp  = r[0] - right_sp;

            return {{{left_sp, left_sp + left_size - 1},
                    {right_sp, right_sp + right_size - 1}
                }
            };
        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const
        {
            return {m_max_level,c};
        }

    private:

        //! Iterator to the begin of the bitvector of inner node v
        auto begin(const node_type& v) const -> decltype(m_tree.begin() + v.offset)
        {
            return m_tree.begin() + v.offset;
        }

        //! Iterator to the begin of the bitvector of inner node v
        auto end(const node_type& v) const -> decltype(m_tree.begin() + v.offset + v.size)
        {
            return m_tree.begin() + v.offset + v.size;
        }
};

}// end namespace sdsl
#endif
