/*! \file cst_fully.hpp
    \brief cst_fully.hpp contains an implementation of Russo et al.'s Fully-Compressed Suffix Tree.
    \author Christian Ocker, Simon Gog
*/
#ifndef INCLUDED_SDSL_CST_FULLY
#define INCLUDED_SDSL_CST_FULLY

#include "bit_vectors.hpp"
#include "bp_support.hpp"
#include "suffix_arrays.hpp"
#include "util.hpp"
#include "vectors.hpp"
#include "cst_sada.hpp"
#include "cst_iterators.hpp"
#include "sdsl_concepts.hpp"
#include "construct.hpp"
#include "suffix_tree_helper.hpp"
#include "suffix_tree_algorithm.hpp"

namespace sdsl
{

template<typename t_cst>
class lcp_fully
{
    public:
        typedef typename t_cst::size_type               size_type;
        typedef size_type                               value_type;
        typedef random_access_const_iterator<lcp_fully> const_iterator;
        typedef const_iterator                          iterator;

        typedef lcp_tag                                 lcp_category;

        enum { fast_access = 0,
               text_order = 0,
               sa_order = 0
             };
    private:
        const t_cst* m_cst;
    public:
        lcp_fully() = default;
        lcp_fully(const t_cst* cst) : m_cst(cst) {};

        lcp_fully(const lcp_fully&) = default;
        lcp_fully(lcp_fully&&) = default;
        lcp_fully& operator=(const lcp_fully&) = default;
        lcp_fully& operator=(lcp_fully&&) = default;
        ~lcp_fully() = default;

        size_type size() const
        {
            return m_cst->size();
        }

        value_type operator[](size_type i) const
        {
            if (0 == i) {
                return 0;
            } else {
                using leaf_type = typename t_cst::leaf_type;
                using char_type = typename t_cst::char_type;
                using sampled_node_type = typename t_cst::sampled_node_type;
                leaf_type v_l = i-1;
                leaf_type v_r = i;

                size_type i;
                sampled_node_type u;
                std::vector<char_type> c(m_cst->delta, 0);
                return m_cst->depth_lca(v_l, v_r, i, u, c);
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
};

//! A class for the Fully-Compressed Suffix Tree (FCST) proposed by Russo et al.
/*!
 * \tparam t_csa        Type of a CSA (member of this type is accessible via
 *                      member `csa`, default class is sdsl::wt).
 * \tparam t_delta      Value of the sampling parameter. Larger values result
 *                      in lower space consumption while requiring more time.
 *                      For `t_delta` = 0, delta = log n log log n is used.
 * \tparam t_s_support  Type of a BPS structure (member accessible via member
 *                      `s_support`, default class is sdsl::bp_support_sada),
 * \tparam t_b          Type of a bit vector for the leaf mapping (member
 *                      accessible via member `b`, default class is
 *                      sdsl::sd_vector),
 * \tparam t_depth      Type of an integer vector for the depth of the sampled
 *                      nodes (member accessible via member `depth_sampling`,
 *                      default class is sdsl::dac_vector),
 * \tparam t_sample_leaves Boolean value indicating whether leaves are to be
 *                         sampled. This is helpful for debugging purposes.
 *
 * It also contains a sdsl::bit_vector which represents the balanced
 * parentheses sequence of the sampled tree. This bit_vector can be accessed
 * via member `s`.
 *
 * A node `v` of the `cst_fully` is represented by an integer `i` which
 * corresponds to the position of the opening parenthesis of the parentheses
 * pair \f$(i,\mu(i))\f$ that corresponds to `v` in `s`.
 *
 * \par Reference
 *  Russo, Lu{\'\i}s and Navarro, Gonzalo and Oliveira, Arlindo L:
 *  Fully Compressed Suffix Trees.
 *  ACM Transactions on Algorithms (TALG), vol. 7, no. 4, p. 53, 2011
 *
 * @ingroup cst
 */
template<class t_csa = csa_wt<>,
         uint32_t t_delta = 0,
         class t_s_support = bp_support_sada<>,
         class t_b = sd_vector<>,
         class t_depth = dac_vector<>,
         bool t_sample_leaves = false
         >
class cst_fully
{
    public:
        typedef cst_dfs_const_forward_iterator<cst_fully> const_iterator;
        typedef typename t_csa::size_type                 size_type;
        typedef t_csa                                     csa_type;
        typedef lcp_fully<cst_fully>                      lcp_type;
        typedef typename t_csa::char_type                 char_type;
        typedef std::pair<size_type, size_type>           node_type; // Nodes are represented by their interval over the CSA
        typedef size_type                                 leaf_type; // Index of a leaf
        typedef size_type                                 sampled_node_type; // Node in the sampled tree represented by its index in s
        typedef t_s_support                               s_support_type;
        typedef t_b                                       b_type;
        typedef typename t_b::select_0_type               b_select_0_type;
        typedef typename t_b::select_1_type               b_select_1_type;
        typedef t_depth                                   depth_type;

        typedef typename t_csa::alphabet_category         alphabet_category;
        typedef cst_tag                                   index_category;

    private:
        size_type        m_delta;
        size_type        m_nodes;
        csa_type         m_csa;
        bit_vector       m_s;
        s_support_type   m_s_support;
        b_type           m_b;
        b_select_0_type  m_b_select0;
        b_select_1_type  m_b_select1;
        depth_type       m_depth;
        lcp_type         m_lcp = lcp_type(this);

        void copy(const cst_fully& cst)
        {
            m_delta     = cst.m_delta;
            m_nodes     = cst.m_nodes;
            m_csa       = cst.m_csa;
            m_s         = cst.m_s;
            m_s_support = cst.m_s_support;
            m_s_support.set_vector(&m_s);
            m_b         = cst.m_b;
            m_b_select0 = cst.m_b_select0;
            m_b_select0.set_vector(&m_b);
            m_b_select1 = cst.m_b_select1;
            m_b_select1.set_vector(&m_b);
            m_depth     = cst.m_depth;
        }

    public:
        const size_type&       delta = m_delta;
        const csa_type&        csa = m_csa;
        const bit_vector&      s = m_s;
        const s_support_type&  s_support = m_s_support;
        const b_type&          b = m_b;
        const b_select_0_type& b_select_0 = m_b_select0;
        const b_select_1_type& b_select_1 = m_b_select1;
        const depth_type&      depth_sampling = m_depth;
        const lcp_type&        lcp = m_lcp;

//! Default constructor
        cst_fully() {}

//! Copy constructor
        cst_fully(const cst_fully& cst)
        {
            copy(cst);
        }

//! Move constructor
        cst_fully(cst_fully&& cst)
        {
            *this = std::move(cst);
        }

//! Construct CST from file_map
        cst_fully(cache_config& config);

        size_type size() const
        {
            return m_csa.size();
        }

        static size_type max_size()
        {
            return t_csa::max_size();
        }

        bool empty() const
        {
            return m_csa.empty();
        }

        void swap(cst_fully& cst)
        {
            if (this != &cst) {
                std::swap(m_delta, cst.m_delta);
                std::swap(m_nodes, cst.m_nodes);
                m_csa.swap(cst.m_csa);
                m_s.swap(cst.m_s);
                util::swap_support(m_s_support, cst.m_s_support, &m_s, &(cst.m_s));
                m_b.swap(cst.m_b);
                util::swap_support(m_b_select0, cst.m_b_select0, &m_b, &(cst.m_b));
                util::swap_support(m_b_select1, cst.m_b_select1, &m_b, &(cst.m_b));
                m_depth.swap(cst.m_depth);
            }
        }

        const_iterator begin() const
        {
            if (m_b.size() == 0) {
                return end();
            }
            return const_iterator(this, root(), false, true);
        }

        const_iterator end() const
        {
            return const_iterator(this, root(), true, false);
        }

//! Copy Assignment Operator.
        cst_fully& operator=(const cst_fully& cst)
        {
            if (this != &cst) {
                copy(cst);
            }
            return *this;
        }

//! Move Assignment Operator.
        cst_fully& operator=(cst_fully &&cst)
        {
            if (this != &cst) {
                m_delta     = cst.m_delta;
                m_nodes     = cst.m_nodes;
                m_csa       = std::move(cst.m_csa);
                m_s         = std::move(cst.m_s);
                m_s_support = std::move(cst.m_s_support);
                m_s_support.set_vector(&m_s);
                m_b         = std::move(cst.m_b);
                m_b_select0 = std::move(cst.m_b_select0);
                m_b_select0.set_vector(&m_b);
                m_b_select1 = std::move(cst.m_b_select1);
                m_b_select1.set_vector(&m_b);
                m_depth     = std::move(cst.m_depth);
            }
            return *this;
        }

//! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="") const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_delta, out, child, "m_delta");
            written_bytes += write_member(m_nodes, out, child, "m_nodes");
            written_bytes += m_csa.serialize(out, child, "csa");
            written_bytes += m_s.serialize(out, child, "s");
            written_bytes += m_s_support.serialize(out, child, "s_support");
            written_bytes += m_b.serialize(out, child, "b");
            written_bytes += m_b_select0.serialize(out, child, "b_select0");
            written_bytes += m_b_select1.serialize(out, child, "b_select1");
            written_bytes += m_depth.serialize(out, child, "depth");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

//! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in)
        {
            read_member(m_delta, in);
            read_member(m_nodes, in);
            m_csa.load(in);
            m_s.load(in);
            m_s_support.load(in, &m_s);
            m_b.load(in);
            m_b_select0.load(in, &m_b);
            m_b_select1.load(in, &m_b);
            m_depth.load(in);
        }

//! Returns the root of the suffix tree.
        node_type root() const
        {
            return node_type(0, m_csa.size() - 1);
        }

//! Returns the root of the sampled tree.
        sampled_node_type sampled_root() const
        {
            return 0;
        }

//! Returns true iff node v is a leaf.
        bool is_leaf(node_type v) const
        {
            return v.first == v.second;
        }

//! Return the i-th leaf (1-based from left to right) of the suffix tree.
        /*!
         * \param i 1-based position of the leaf. \f$1\leq i\leq csa.size()\f$.
         * \return The i-th leave.
         * \par Time complexity
         *     \f$ \Order{1} \f$
         * \pre \f$ 1 \leq i \leq csa.size() \f$
         */
        node_type select_leaf(size_type i) const
        {
            assert(i > 0 and i <= m_csa.size());
            return node_type(i - 1, i - 1);
        }

//! Get the node in the suffix tree which corresponds to the sa-interval [lb..rb]
        node_type node(size_type lb, size_type rb) const
        {
            return node_type(lb, rb);
        }

//! Returns the leftmost leaf (left boundary) of a node.
        leaf_type lb(node_type v) const
        {
            return v.first;
        }

//! Returns the rightmost leaf (right boundary) of a node.
        leaf_type rb(node_type v) const
        {
            return v.second;
        }

        //! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The number of leaves in the subtree rooted at node v.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        size_type size(const node_type& v) const
        {
            return v.second-v.first+1;
        }

//! Calculates the leftmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The leftmost leaf in the subtree rooted at node v.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        node_type leftmost_leaf(const node_type v) const
        {
            return node_type(v.first, v.first);
        }

//! Calculates the rightmost leaf in the subtree rooted at node v.
        /*!\param v A valid node of the suffix tree.
         * \return The rightmost leaf in the subtree rooted at node v.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        node_type rightmost_leaf(const node_type v) const
        {
            return node_type(v.second, v.second);
        }

//! Returns true iff v is an ancestor of w.
        bool ancestor(node_type v, node_type w) const
        {
            return v.first <= w.first && v.second >= w.second;
        }

//! Returns the index of the last bracket in S before the leaf with index l.
        /*!
         * \param v The index of leaf l.
         * \return The index of the last bracket in S before the leaf with index l.
         */
        sampled_node_type pred(leaf_type v) const
        {
            return m_b_select0.select(v + 1) - v - 1;
        }

//! Returns the LSA (lowest sampled ancestor) for a leaf with index l.
        /*!
         * \param v The index of leaf l.
         * \return The LSA for the leaf with index l.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        sampled_node_type lsa_leaf(leaf_type l) const
        {
            sampled_node_type p = pred(l);
            if (m_s[p]) {
                return p;
            } else {
                return m_s_support.enclose(m_s_support.find_open(p));
            }
        }

//! Returns the node in the suffix tree corresponding to the node u in the sampled tree.
        /*!
         * \param v The node u in the sampled tree.
         * \return The node in the suffix tree corresponding to the node u in the sampled tree.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        node_type sampled_node(sampled_node_type u) const
        {
            assert(m_s[u] == 1);
            size_type u_end = m_s_support.find_close(u);
            size_type b_left = m_b_select1.select(u + 1);
            size_type b_right = m_b_select1.select(u_end + 1);

            return node_type(b_left - u,
                             b_right - u_end - 1);
        }

//! Returns the LCA of two nodes in the sampled tree.
        /*!
         * \param u The sampled node u.
         * \param q The sampled node q.
         * \return The lowest common ancestor of u and q in the sampled tree.
         * \par Time complexity
         *   \f$ \Order{\rrenclose} \f$
         */
        sampled_node_type sampled_lca(sampled_node_type u, sampled_node_type q) const
        {
            assert(m_s[u] == 1 and m_s[q] == 1);
            if (u > q) {
                std::swap(u, q);
            } else if (u == q) {
                return u;
            }
            if (u == sampled_root()) {
                return sampled_root();
            }
            if (m_s_support.find_close(u) > q) {
                return u;
            }

            return m_s_support.double_enclose(u, q);
        }

//! Returns the depth of a sampled node u.
        /*!
         * \param u A sampled node u.
         * \return The depth of sampled node u.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        size_type depth(sampled_node_type u) const
        {
            assert(m_s[u] == 1);

            size_type idx = m_s_support.rank(u) - 1;
            return m_depth[idx] * (m_delta / 2);
        }

//! Returns the depth of a node v.
        /*!
         * \param v The node v.
         * \return The depth of node v.
         * \par Time complexity
         *   \f$ \Order( \delta ) \f$ for inner nodes,
         *   \f$ \Order( \saaccess ) \f$ for leaves.
         */
        size_type depth(node_type v) const
        {
            if (is_leaf(v)) {
                return m_csa.size() - m_csa[v.first];
            }

            size_type i;
            sampled_node_type u;
            std::vector<char_type> c;
            c.reserve(delta);
            return depth_lca(v.first, v.second, i, u, c);
        }

//! Calculate the LCA of two nodes v and w.
        /*!
         * \param v The node v.
         * \param w The node w.
         * \return The LCA of v and w.
         * \par Time complexity
         *   \f$ \Order( \delta \cdot ( 1 + t_{rank\_bwt} ) ) \f$
         */
        node_type lca(node_type v, node_type w) const
        {
            leaf_type l = std::min(v.first, w.first);
            leaf_type r = std::max(v.second, w.second);

            if (l == r) {
                return node_type(l, r);
            } else {
                return lca(l, r);
            }
        }

//! Calculate the LCA of two leaves l and r.
        /*!
         * \param l The index of leaf l.
         * \param r The index of leaf r. \f$ r > l \f$
         * \return The LCA of l and r.
         * \par Time complexity
         *   \f$ \Order( \delta \cdot ( 1 + t_{rank\_bwt} ) ) \f$
         */
        node_type lca(leaf_type l, leaf_type r) const
        {
            assert(l<r);

            size_type i;
            sampled_node_type u;
            std::vector<char_type> c(delta, 0);
            depth_lca(l, r, i, u, c);

            node_type v = sampled_node(u);
            leaf_type lb = v.first;
            leaf_type rb = v.second;

            for (size_type k = 0; k < i; k++) {
                backward_search(m_csa, lb, rb, c[i - k - 1], lb, rb);
            }

            return node_type(lb, rb);
        }

//! Calculate the depth of the LCA of two leaves l and r.
        /*!
         * \param l The index of leaf l.
         * \param r The index of leaf r. \f$ r > l \f$
         * \param res_i The index i for the ancestor used to determine the depth (return value).
         * \param res_u The ancestor used to determine the depth (return value).
         * \param res_label The label from the found sampled node to the actual LCA.
         * \return The depth of the LCA of l and r.
         * \par Time complexity
         *   \f$ \Order( \delta ) \f$
         */
        // TODO: return by reference really necessary?
        size_type depth_lca(leaf_type l, leaf_type r,
                            size_type& res_i, sampled_node_type& res_u, std::vector<char_type>& res_label) const
        {
            assert(l<r);

            size_type max_d = 0;
            size_type max_d_i = 0;
            sampled_node_type max_d_node = 0;

            for (size_type i = 0; i < m_delta; i++) {
                sampled_node_type node = sampled_lca(lsa_leaf(l), lsa_leaf(r));
                size_type d = i + depth(node);

                if (d > max_d) {
                    max_d = d;
                    max_d_i = i;
                    max_d_node = node;
                }

                char_type c = m_csa.F[l];
                char_type comp = csa.char2comp[c];
                res_label[i] = c;

                // break if LCA of lb and rb is root
                if (l < m_csa.C[comp] || r >= m_csa.C[comp + 1]) {
                    break;
                }

                l = m_csa.psi[l];
                r = m_csa.psi[r];
            }

            res_i = max_d_i;
            res_u = max_d_node;

            return max_d;
        }

//! Compute the suffix link of a node v.
        /*!
         * \param v The node v.
         * \return The suffix link of node v or root() if v equals root().
         * \par Time complexity
         *   \f$ \Order( \delta \cdot ( 1 + t_{rank\_bwt} ) ) \f$
         */
        node_type sl(node_type v) const
        {
            if (v == root()) {
                return root();
            } else if (is_leaf(v)) {
                size_t leaf = m_csa.psi[v.first];
                return node_type(leaf, leaf);
            }

            return lca(m_csa.psi[v.first], m_csa.psi[v.second]);
        }

//! Compute the Weiner link of node v and character c.
        /*
         * \param v A valid node of a cst_fully.
         * \param c The character which should be prepended to the string of the current node.
         *   \return root() if the Weiner link of (v, c) does not exist, otherwise the Weiner link is returned.
         * \par Time complexity
         *    \f$ \Order{ t_{rank\_bwt} + t_{lca}}\f$
         */
        node_type wl(node_type v, const char_type c) const
        {
            size_type l, r;
            std::tie(l, r) = v;
            backward_search(m_csa, l, r, c, l, r);
            return node_type(l, r);
        }

//! Compute the suffix number of a leaf node v.
        /*!\param v A valid leaf node of a cst_sada.
         * \return The suffix array value corresponding to the leaf node v.
         * \par Time complexity
         *   \f$ \Order{ \saaccess } \f$
         */
        size_type sn(node_type v) const
        {
            assert(is_leaf(v));
            return m_csa[v.first];
        }

//! Calculate the parent node of a node v.
        /*!
         * \param v The node v.
         * \return The parent node of v or root() if v equals root().
         * \par Time complexity
         *   \f$ \Order( \delta \cdot ( 1 + t_{rank\_bwt} ) ) \f$
         */
        node_type parent(node_type v) const
        {
            const leaf_type l = v.first;
            const leaf_type r = v.second;

            node_type left_parent = root();
            node_type right_parent = root();

            if (l > 0) {
                left_parent = lca(l-1, r);
            }
            if (r < m_csa.size() - 1) {
                right_parent = lca(l, r+1);
            }
            return ancestor(right_parent, left_parent) ? left_parent : right_parent;
        }

//! Get the child w of node v which edge label (v,w) starts with character c.
        /*!
         * \param v A node v.
         * \param c First character of the edge label from v to the desired child.
         * \return The child node w which edge label (v,w) starts with c or root() if it does not exist.
         * \par Time complexity
         *       \f$ \Order{ \log m \cdot (\saaccess+\isaaccess) } \f$
                 where \f$ m \f$ is the number of leaves in the subtree rooted at node v.
         */
        node_type child(node_type v, char_type c) const
        {
            if (is_leaf(v)) {
                return root();
            }
            size_type d = depth(v);
            return child(v, c, d);
        }

        node_type child(node_type v, char_type c, size_type d) const
        {
            leaf_type lower;
            leaf_type upper;

            {
                leaf_type begin = v.first;
                leaf_type end = v.second + 1;

                while (begin < end) {
                    leaf_type sample_pos = (begin + end) / 2;
                    size_type char_pos = get_char_pos(sample_pos, d, m_csa);
                    char_type sample = m_csa.F[char_pos];
                    if (sample < c) {
                        begin = sample_pos + 1;
                    } else {
                        end = sample_pos;
                    }
                }

                lower = begin;
            }

            {
                leaf_type begin = v.first;
                leaf_type end = v.second + 1;

                while (begin < end) {
                    leaf_type sample_pos = (begin + end) / 2;
                    size_type char_pos = get_char_pos(sample_pos, d, m_csa);
                    char_type sample = m_csa.F[char_pos];
                    if (sample <= c) {
                        begin = sample_pos + 1;
                    } else {
                        end = sample_pos;
                    }
                }

                upper = begin;
            }

            if (lower == upper) {
                return root();
            }

            return node_type(lower, upper - 1);
        }

//! Get the i-th child of a node v.
        /*!
         * \param v A valid tree node of the cst.
         * \param i 1-based Index of the child which should be returned. \f $i \geq 1 \f$.
         * \return The i-th child node of v or root() if v has no i-th child.
         */
        node_type select_child(node_type v, size_type i) const
        {
            if (is_leaf(v)) {
                return root();
            }

            size_type d = depth(v);
            size_type char_pos = get_char_pos(v.first, d, m_csa);
            char_type c = m_csa.F[char_pos];
            node_type res = child(v, c, d);
            while (i > 1) {
                if (res.second >= v.second) {
                    return root();
                }
                char_pos = get_char_pos(res.second + 1, d, m_csa);
                c = m_csa.F[char_pos];
                res = child(v, c, d);
                i--;
            }

            return res;
        }

        //! Get the number of children of a node v.
        /*!
         * \param v A valid node v.
         * \returns The number of children of node v.
         */
        size_type degree(const node_type& v)const
        {
            if (is_leaf(v)) {
                return 0;
            } else {
                size_type res = 1;
                size_type d = depth(v);
                size_type char_pos = get_char_pos(v.first, d, m_csa);
                char_type c = m_csa.F[char_pos];
                node_type v_i = child(v, c, d);
                while (v_i.second < v.second) {
                    ++res;
                    char_pos = get_char_pos(v_i.second + 1, d, m_csa);
                    c = m_csa.F[char_pos];
                    v_i = child(v, c, d);
                }
                return res;
            }
        }

        //! Return a proxy object which allows iterating over the children of a node
        /*! \param v A valid node of the suffix tree.
         *  \return The proxy object of v containing all children
         */
        cst_node_child_proxy<cst_fully> children(const node_type& v) const
        {
            return cst_node_child_proxy<cst_fully>(this,v);
        }



//! Returns the next sibling of node v.
        /*!
         * \param v A valid node v of the suffix tree.
         * \return The next (right) sibling of node v or root() if v has no next sibling.
         */
        node_type sibling(node_type v) const
        {
            node_type p = parent(v);
            if (v.second >= p.second) {
                return root();
            }
            size_type d = depth(p);
            size_type char_pos = get_char_pos(v.second + 1, d, m_csa);
            char_type c = m_csa.F[char_pos];
            return child(p, c, d);
        }

        char_type edge(node_type v, size_type d) const
        {
            assert(d >= 1 and d <= depth(v));
            size_type char_pos = get_char_pos(v.first, d - 1, m_csa);
            return m_csa.F[char_pos];
        }

        //! Returns the node depth of node v
        /*!
         * \param v A valid node of a cst_fully
         * \return The node depth of node v.
         */
        size_type node_depth(node_type v)const
        {
            size_type d = 0;
            while (v != root()) {
                ++d;
                v = parent(v);
            }
            return d;
        }

        //! Get the number of nodes of the suffix tree.
        size_type nodes()const
        {
            return m_nodes;
        }

//! Get the number of nodes in the sampled tree.
        /*!
         * \return The number of nodes in the sampled tree.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        size_type sampled_nodes() const
        {
            return m_s.size() / 2;
        }
};

template<class t_csa, uint32_t t_delta, class t_s_support, class t_b, class t_depth, bool t_sample_leaves>
cst_fully<t_csa, t_delta, t_s_support, t_b, t_depth, t_sample_leaves>::cst_fully(cache_config& config)
{
    // 1. Construct CST
    cst_sada<csa_type, lcp_dac<> > cst(config);
    m_nodes = cst.nodes();

    if (t_delta > 0) {
        m_delta = t_delta;
    } else {
        const size_type n = cst.size();
        m_delta = (bits::hi(n-1)+1) * (bits::hi(bits::hi(n-1))+1);
        if (m_delta < 2) {
            m_delta = 2;
        }
    }

    size_type delta_half = m_delta / 2;

    bit_vector is_sampled(cst.nodes(), false);
    is_sampled[cst.id(cst.root())] = true; // always sample root
    size_type sample_count = 1;

    // 2a. Scan and mark leaves to be sampled
    if (t_sample_leaves) {
        auto event = memory_monitor::event("scan-leaves");
        size_type leaf_idx = 0;
        for (size_type i = 0; i < cst.size(); i++) {
            const size_type d = i + 1;
            if (d + delta_half <= cst.size() and
                d % delta_half == 0) {
                const auto node = cst.select_leaf(leaf_idx + 1);
                const size_type id = cst.id(node);
                if (!is_sampled[id]) {
                    is_sampled[id] = true;
                    sample_count++;
                }
            }
            leaf_idx = cst.csa.lf[leaf_idx];
        }
    }

    // 2b. Scan and mark inner nodes to be sampled
    {
        auto event = memory_monitor::event("scan-nodes");
        for (auto it = cst.begin(); it != cst.end(); ++it) {
            if (it.visit() == 1 and cst.is_leaf(*it) == false) {
                const auto node = *it;
                const size_type d = cst.depth(node);
                if (d % delta_half == 0) {
                    auto v = cst.sl(node, delta_half);
                    const size_type id = cst.id(v);
                    if (!is_sampled[id]) {
                        is_sampled[id] = true;
                        sample_count++;
                    }
                }
            }
        }
    }

    m_s.resize(2 * sample_count);
    util::set_to_value(m_s, 0);
    bit_vector tmp_b(2 * sample_count + cst.size(), 0);
    int_vector<64> tmp_depth;
    tmp_depth.resize(sample_count);

    // 3. Create sampled tree data structures
    {
        auto event = memory_monitor::event("node-sampling");

        size_type s_idx = 0;
        size_type b_idx = 0;
        size_type sample_idx = 0;

        for (auto it = cst.begin(); it != cst.end(); ++it) {
            auto node = *it;
            if (it.visit() == 1 && is_sampled[cst.id(node)]) {
                m_s[s_idx++] = 1;
                tmp_b[b_idx++] = 1;
                tmp_depth[sample_idx++] = cst.depth(node) / delta_half;
            }
            if (cst.is_leaf(node)) {
                b_idx++;
            }
            if ((cst.is_leaf(node) || it.visit() == 2) && is_sampled[cst.id(node)]) {
                s_idx++;
                tmp_b[b_idx++] = 1;
            }
        }
    }

    {
        auto event = memory_monitor::event("ss-depth");
        m_csa = std::move(cst.csa);
        util::init_support(m_s_support, &m_s);
        m_b = b_type(tmp_b);
        util::init_support(m_b_select0, &m_b);
        util::init_support(m_b_select1, &m_b);
        m_depth = depth_type(tmp_depth);
    }


}

}// end namespace sdsl

// TODO: make dependent on cst_fully
template<class T>
std::ostream& operator<<(std::ostream& os, const std::pair<T, T>& v)
{
    os << "[" << v.first << ", " << v.second << "]";
    return os;
}

#endif // INCLUDED_SDSL_CST_FULLY
