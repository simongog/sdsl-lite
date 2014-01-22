#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;
using namespace std::chrono;


template<class t_bitvector = bit_vector,
         class t_select = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
class wt_gmr_1
{
    private:

        t_bitvector m_bv;
        int_vector<> e;
        t_select sls1;
        t_select_zero sls0;
        uint64_t m_size; // input length
        uint64_t m_sigma = 0; // maximum character + 1
        uint64_t m_blocks; // blocks per character

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
                if (m_sigma < input[i]) m_sigma = input[i];
            }
            ++m_sigma;

            // Create and fill b
            m_blocks = (m_size+m_sigma-1)/m_sigma;
            bit_vector b(m_size+m_sigma*m_blocks);
            int_vector<> symbols(m_sigma,0,bits::hi(m_size)+1);
            {
                int_vector<> tmp(m_sigma*m_blocks,0,bits::hi(m_sigma)+1);

                for (uint64_t i=0, offset=0, j=0; i<m_size; ++i, ++j) {
                    if (j==m_sigma) {
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
                /*
                for(uint64_t i=0, l=0; i<tmp.size(); ++i, ++l){
                	while(tmp[i] > 64) {
                		b.set_int(l, 0xFFFFFFFFFFFFFFFFULL, 64);
                		l += 64;
                		tmp[i] -= 64;
                	}
                	b.set_int(l, 0xFFFFFFFFFFFFFFFFULL, tmp[i]);
                	l += tmp[i];
                }
                */
                for (uint64_t i=0,l=0; i<tmp.size(); ++i,++l) {
                    for (uint64_t j=0; j<tmp[i]; ++j) b[l++]=1;
                }
                m_bv = t_bitvector(std::move(b));
            }
            util::init_support(sls0, &m_bv);
            util::init_support(sls1, &m_bv);

            // Create and fill e
            e = int_vector<>(m_size,0,bits::hi(m_sigma)+1);
            for (uint64_t i=0, tmp=0, sum=0; i<m_sigma; ++i) {
                tmp = symbols[i];
                symbols[i] = sum;
                sum += tmp;
            }
            for (uint64_t i=0; i<m_size;) {
                for (uint64_t j=0; j<m_sigma and i<m_size; ++i, ++j) {
                    e[symbols[input[i]]++] = j;
                }
            }
        }

        //! Swap operator
        void swap(wt_gmr_1& fs) {
            if (this != &fs) {
                m_bv.swap(fs.m_bv);
                e.swap(fs.e);
                util::swap_support(sls0, fs.sls0, &m_bv, &(fs.m_bv));
                util::swap_support(sls1, fs.sls1, &m_bv, &(fs.m_bv));
                std::swap(m_size, fs.m_size);
                std::swap(m_sigma,  fs.m_sigma);
                std::swap(m_blocks,  fs.m_blocks);
            }
        }

        value_type operator[](size_type i)const {
            return i;
        }

        pair<size_type, value_type>	inverse_select(size_type i)const {
            return make_pair(i,0);
        }

        size_type select(size_type i, value_type c)const {
            size_type k = (c==0?0:sls0(c*m_blocks)-c*m_blocks+1)+i;
            return (sls1(k)-k+1)*m_sigma+e[k-1]-c*m_blocks*m_sigma;
        }

        size_type rank(size_type i, value_type c)const {
            if (c>m_sigma-1) return 0;
            if (i<=0) return 0;

            size_type ones_before_cblock;
            size_type search_begin;
            size_type search_end;
            size_type offset=0;
            if (c!=0) {
                ones_before_cblock = sls0(c*m_blocks)-c*m_blocks+1;
                search_begin = sls0(c*m_blocks+(i-1)/m_sigma)-(c*m_blocks+(i-1)/m_sigma)+1;
                search_end = sls0(c*m_blocks+(i-1)/m_sigma+1)-(c*m_blocks+(i-1)/m_sigma+1)+1;
            } else {
                ones_before_cblock = 0;
                search_begin = ((i-1)<m_sigma?0:sls0((i-1)/m_sigma)-(i-1)/m_sigma+1);
                search_end = sls0((i-1)/m_sigma+1)-((i-1)/m_sigma+1)+1;
                if (search_end-search_begin==0) return 0;
            }

            size_type val = (i-1)%m_sigma;
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
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += write_member(m_blocks, out, child, "blocks");
            written_bytes += m_bv.serialize(out, child, "b");
            written_bytes += e.serialize(out, child, "e");
            written_bytes += sls0.serialize(out, child, "sls0");
            written_bytes += sls1.serialize(out, child, "sls1");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

template<class t_bitvector = bit_vector,
         class t_select = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_rank = typename t_bitvector::rank_1_type>
class wt_gmr_2
{
    private:

        t_bitvector m_bv;
        t_bitvector m_xv;
        t_bitvector m_piv;

        int_vector<> s;
        int_vector<> pi;

        t_rank piv_r1;
        t_select b_sls1, x_sls1;
        t_select_zero b_sls0, x_sls0;

        uint64_t m_size; // input length
        uint64_t m_t = 32; // shortcut value
        uint64_t m_sigma = 0; // maximum character + 1
        uint64_t m_chunks; // number of chunks

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
                if (m_sigma < input[i]) m_sigma = input[i];
            }
            ++m_sigma;

            m_chunks = (m_size+m_sigma-1)/m_sigma;


            pi = int_vector<>(m_size,0,bits::hi(m_sigma-1)+1);
            {
                uint64_t x_pos = 0;
                bit_vector x(m_sigma*m_chunks+m_size);

                //fill pi and x for every chunk
                for (uint64_t i=0; i<m_chunks; ++i) {
                    int_vector<> symbols(m_sigma,0,bits::hi(m_sigma-1)+1);

                    //calc symbols
                    for (uint64_t j=i*m_sigma; j<(i+1)*m_sigma and j<m_size; ++j) {
                        ++symbols[input[j]];
                    }
                    //calc x
                    for (uint64_t j=0; j<m_sigma; ++j,++x_pos) {
                        for (uint64_t k=0; k<symbols[j]; ++k) x[++x_pos]=1;
                    }

                    //calc symbols prefix sum
                    for (uint64_t j=0, tmp=0, sum=0; j<m_sigma; ++j) {
                        tmp = symbols[j];
                        symbols[j] = sum;
                        sum += tmp;
                    }
                    //calc pi
                    for (uint64_t j=i*m_sigma, k=0; j<(i+1)*m_sigma and j<m_size; ++j,++k) {
                        pi[i*m_sigma+(symbols[input[j]]++)] = k;
                    }
                }
                m_xv = t_bitvector(std::move(x));
                util::init_support(x_sls1, &m_xv);
                util::init_support(x_sls0, &m_xv);

            }
            //calc b
            {
                bit_vector b(m_size+m_sigma*m_chunks,0);
                int_vector<> tmp(m_sigma*m_chunks,0,bits::hi(m_sigma-1)+1);

                for (uint64_t i=0, offset=0, j=0; i<m_size; ++i, ++j) {
                    if (j==m_sigma) {
                        ++offset;
                        j = 0;
                    }
                    ++tmp[offset+input[i]*m_chunks];
                }

                for (uint64_t i=0,l=0; i<tmp.size(); ++i,++l) {
                    for (uint64_t j=0; j<tmp[i]; ++j) b[l++]=1;
                }
                m_bv = t_bitvector(std::move(b));
                util::init_support(b_sls1, &m_bv);
                util::init_support(b_sls0, &m_bv);
            }
            //calc inverse pi
            {
                bit_vector ipi(m_size);
                bit_vector pitmp(m_sigma,0);
                //calc pointer pos
                for (uint64_t i=0; i<m_chunks; ++i) {
                    bit_vector pitmp(m_sigma,0);
                    for (uint64_t j = i*m_sigma,k=0; j<(i+1)*m_sigma and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_sigma,cycle_length=0;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    ipi[pos_pi]=1;
                                    steps = 0;
                                }
                                pos_pitmp = pi[pos_pi];
                                pos_pi = offset + pi[pos_pi];
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
                m_piv = t_bitvector(std::move(ipi));
                util::init_support(piv_r1, &m_piv);

                //calc s
                s= int_vector<>(piv_r1(m_size),0,bits::hi(m_sigma-1)+1);

                for (uint64_t i=0; i<m_chunks; ++i) {
                    //bit_vector pitmp(m_sigma,0);
                    for (uint64_t j = i*m_sigma,k=0; j<(i+1)*m_sigma and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_sigma,cycle_length=0,back_pointer=k;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    s[piv_r1(pos_pi)]=back_pointer;
                                    back_pointer=pos_pitmp;
                                    steps = 0;
                                }
                                pos_pitmp = pi[pos_pi];
                                pos_pi = offset + pi[pos_pi];
                                ++steps;
                                ++cycle_length;
                            } while (!pitmp[pos_pitmp]);

                            if (steps==m_t and cycle_length>m_t) {
                                s[piv_r1(pos_pi)]=back_pointer;
                            }
                        }
                    }
                    util::_set_zero_bits(pitmp);
                }
            }
        }



        //add members
        //! Swap operator
        void swap(wt_gmr_2& fs) {
            if (this != &fs) {
                m_bv.swap(fs.m_bv);
                m_xv.swap(fs.m_xv);
                m_piv.swap(fs.m_piv);
                pi.swap(fs.pi);
                s.swap(fs.s);
                util::swap_support(b_sls0, fs.b_sls0, &m_bv, &(fs.m_bv));
                util::swap_support(b_sls1, fs.b_sls1, &m_bv, &(fs.m_bv));
                util::swap_support(x_sls1, fs.x_sls1, &m_xv, &(fs.m_xv));
                util::swap_support(x_sls0, fs.x_sls0, &m_xv, &(fs.m_xv));
                util::swap_support(piv_r1, fs.piv_r1, &m_piv, &(fs.m_piv));
                std::swap(m_size, fs.m_size);
                std::swap(m_t, fs.m_t);
                std::swap(m_sigma,  fs.m_sigma);
                std::swap(m_chunks,  fs.m_chunks);
            }
        }

        value_type operator[](size_type i)const {
            uint64_t chunk = i/m_sigma;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_sigma);
            while (pi[x]!=value) {
                if (jump and m_piv[x]==1) {
                    x  = s[piv_r1(x)]+(chunk*m_sigma);
                    jump = false;
                } else {
                    x = pi[x]+(chunk*m_sigma);
                }
            }
            return x_sls1(x+1)-x-(chunk*m_sigma)-1;
        }

        pair<size_type, value_type>	inverse_select(size_type i)const {
            uint64_t chunk = i/m_sigma;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_sigma);
            while (pi[x]!=value) {
                if (jump and m_piv[x]==1) {
                    x  = s[piv_r1(x)]+(chunk*m_sigma);
                    jump = false;
                } else {
                    x = pi[x]+(chunk*m_sigma);
                }
            }
            uint64_t tmp = x_sls1(x+1);
            uint64_t c = tmp-x-(chunk*m_sigma)-1;
            uint64_t c_before_chunk;
            if (c==0) {
                if (chunk==0) {
                    c_before_chunk =0;
                } else {
                    c_before_chunk = b_sls0(chunk)-chunk+1;
                }
            } else {
                uint64_t ones_before_c = b_sls0(c*m_chunks)-(c*m_chunks)+1;
                c_before_chunk = b_sls0(c*m_chunks+chunk)-(c*m_chunks+chunk)+1-ones_before_c;
            }
            uint64_t c_in_chunk = tmp-x_sls0(c+1+chunk*m_sigma);
            return make_pair(c_before_chunk+c_in_chunk,c);
        }

        size_type select(size_type i, value_type c)const {

            uint64_t chunk, c_ones_before_chunk, pi_pos;

            if (c==0) {
                chunk = b_sls1(i)-i+1;
                if (chunk==0) {
                    c_ones_before_chunk = 0;
                } else {
                    c_ones_before_chunk = b_sls0(chunk)-chunk+1;
                }
                pi_pos = i-c_ones_before_chunk-1+chunk*m_sigma;
            } else {
                uint64_t ones_before_c = b_sls0(c*m_chunks)-(c*m_chunks)+1;
                chunk = b_sls1(ones_before_c+i)-ones_before_c-(c*m_chunks)-i+1;
                c_ones_before_chunk = b_sls0(c*m_chunks+chunk)-(c*m_chunks+chunk)+1-ones_before_c;
                pi_pos = x_sls0(chunk*m_sigma+c+1)+(i-c_ones_before_chunk)-chunk*m_sigma-c-1;
            }

            return pi[pi_pos]+chunk*m_sigma;
        }

        size_type rank(size_type i, value_type c)const {

            if (c>m_sigma-1) return 0;
            if (i<=0) return 0;

            uint64_t chunk = (i-1)/m_sigma;
            uint64_t c_ones_before_chunk;

            if (c==0) {
                if (chunk==0) {
                    c_ones_before_chunk=0;
                } else {
                    c_ones_before_chunk=b_sls0(chunk)-chunk+1;
                }
            } else {
                uint64_t ones_before_c = b_sls0(c*m_chunks)-(c*m_chunks)+1;
                c_ones_before_chunk = b_sls0(c*m_chunks+chunk)-(c*m_chunks+chunk)+1-ones_before_c;
            }
            uint64_t c_ones_in_chunk = 0;

            size_type search_begin = x_sls0(chunk*m_sigma+1+c)-(chunk*m_sigma+1+c)+1;
            size_type search_end = x_sls0(chunk*m_sigma+2+c)-(chunk*m_sigma+2+c)+1;

            size_type val = (i-1)%m_sigma;
            if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                while (search_begin < search_end and pi[search_begin] <= val) {
                    ++search_begin;
                    ++c_ones_in_chunk;
                }
            } else {
                c_ones_in_chunk = lower_bound(pi.begin()+search_begin, pi.begin()+search_end, val+1)-pi.begin()-search_begin;
            }
            return c_ones_before_chunk+c_ones_in_chunk;
        }

        //add members
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += write_member(m_chunks, out, child, "chunks");
            written_bytes += write_member(m_t, out, child, "t");
            written_bytes += m_bv.serialize(out, child, "b");
            written_bytes += m_xv.serialize(out, child, "x");
            written_bytes += m_piv.serialize(out, child, "piv");
            written_bytes += pi.serialize(out, child, "pi");
            written_bytes += s.serialize(out, child, "s");
            written_bytes += b_sls0.serialize(out, child, "b_sls0");
            written_bytes += b_sls1.serialize(out, child, "b_sls1");
            written_bytes += x_sls1.serialize(out, child, "x_sls0");
            written_bytes += x_sls0.serialize(out, child, "x_sls1");
            written_bytes += piv_r1.serialize(out, child, "piv_r1");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

//tmp
template<class t_bitvector = bit_vector,
         class t_select = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_rank = typename t_bitvector::rank_1_type>
class wt_gmr_3
{
    private:

        t_bitvector m_bv;
        t_bitvector m_xv;
        t_bitvector m_piv;

        int_vector<> s;
        int_vector<> pi;

        t_rank piv_r1;
        t_select b_sls1, x_sls1;
        t_select_zero b_sls0, x_sls0;

        uint64_t m_size; // input length
        uint64_t m_t = 32; // shortcut value
        uint64_t m_sigma = 0; // maximum character + 1
        uint64_t m_chunks; // number of chunks
        uint64_t m_chunksize;

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
                if (m_sigma < input[i]) m_sigma = input[i];
            }
            ++m_sigma;

            m_chunksize = (1 << (bits::hi(m_sigma-1)+1));
            m_chunks = (m_size+m_chunksize-1)/m_chunksize;

            pi = int_vector<>(m_size,0,bits::hi(m_sigma-1)+1);
            {
                uint64_t x_pos = 0;
                bit_vector x(m_size+m_chunks*m_sigma,0);

                //fill pi and x for every chunk
                for (uint64_t i=0; i<m_chunks; ++i) {
                    int_vector<> symbols(m_sigma,0,bits::hi(m_sigma-1)+2);

                    //calc symbols
                    for (uint64_t j=i*m_chunksize; j<(i+1)*m_chunksize and j<m_size; ++j) {
                        ++symbols[input[j]];
                    }
                    //calc x
                    for (uint64_t j=0; j<m_sigma; ++j,++x_pos) {
                        for (uint64_t k=0; k<symbols[j]; ++k) x[++x_pos]=1;
                    }

                    //calc symbols prefix sum
                    for (uint64_t j=0, tmp=0, sum=0; j<m_sigma; ++j) {
                        tmp = symbols[j];
                        symbols[j] = sum;
                        sum += tmp;
                    }
                    //calc pi
                    for (uint64_t j=i*m_chunksize, k=0; j<(i+1)*m_chunksize and j<m_size; ++j,++k) {
                        pi[i*m_chunksize+(symbols[input[j]]++)] = k;
                    }
                }
                m_xv = t_bitvector(std::move(x));
                util::init_support(x_sls1, &m_xv);
                util::init_support(x_sls0, &m_xv);

            }
            //calc b
            {
                bit_vector b(m_size+m_sigma*m_chunks,0);
                int_vector<> tmp(m_sigma*m_chunks,0,bits::hi(m_sigma-1)+2);

                for (uint64_t i=0, offset=0, j=0; i<m_size; ++i, ++j) {
                    if (j==m_chunksize) {
                        ++offset;
                        j = 0;
                    }
                    ++tmp[offset+input[i]*m_chunks];
                }

                for (uint64_t i=0,l=0; i<tmp.size(); ++i,++l) {
                    for (uint64_t j=0; j<tmp[i]; ++j) b[l++]=1;
                }
                m_bv = t_bitvector(std::move(b));
                util::init_support(b_sls1, &m_bv);
                util::init_support(b_sls0, &m_bv);
            }

            //calc inverse pi
            {
                bit_vector ipi(m_size,0);
                bit_vector pitmp(m_chunksize,0);
                //calc pointer pos
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
                                pos_pitmp = pi[pos_pi];
                                pos_pi = offset + pi[pos_pi];
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
                m_piv = t_bitvector(std::move(ipi));
                util::init_support(piv_r1, &m_piv);

                //calc s
                s= int_vector<>(piv_r1(m_size),0,bits::hi(m_sigma-1)+1);

                for (uint64_t i=0; i<m_chunks; ++i) {
                    for (uint64_t j = i*m_chunksize,k=0; j<(i+1)*m_chunksize and j<m_size; ++j,++k) {
                        if (!pitmp[k]) {
                            uint64_t steps =0,pos_pi=j,pos_pitmp=k,offset =i*m_chunksize,cycle_length=0,back_pointer=k;
                            do {
                                pitmp[pos_pitmp]=1;
                                if (steps==m_t) {
                                    s[piv_r1(pos_pi)]=back_pointer;
                                    back_pointer=pos_pitmp;
                                    steps = 0;
                                }
                                pos_pitmp = pi[pos_pi];
                                pos_pi = offset + pi[pos_pi];
                                ++steps;
                                ++cycle_length;
                            } while (!pitmp[pos_pitmp]);

                            if (steps==m_t and cycle_length>m_t) {
                                s[piv_r1(pos_pi)]=back_pointer;
                            }
                        }
                    }
                    util::_set_zero_bits(pitmp);
                }
            }
        }



        //add members
        //! Swap operator
        void swap(wt_gmr_3& fs) {
            if (this != &fs) {
                m_bv.swap(fs.m_bv);
                m_xv.swap(fs.m_xv);
                m_piv.swap(fs.m_piv);
                pi.swap(fs.pi);
                s.swap(fs.s);
                util::swap_support(b_sls0, fs.b_sls0, &m_bv, &(fs.m_bv));
                util::swap_support(b_sls1, fs.b_sls1, &m_bv, &(fs.m_bv));
                util::swap_support(x_sls1, fs.x_sls1, &m_xv, &(fs.m_xv));
                util::swap_support(x_sls0, fs.x_sls0, &m_xv, &(fs.m_xv));
                util::swap_support(piv_r1, fs.piv_r1, &m_piv, &(fs.m_piv));
                std::swap(m_size, fs.m_size);
                std::swap(m_t, fs.m_t);
                std::swap(m_sigma,  fs.m_sigma);
                std::swap(m_chunks,  fs.m_chunks);
                std::swap(m_chunksize,  fs.m_chunksize);
            }
        }

        value_type operator[](size_type i)const {
            uint64_t chunk = i/m_chunksize;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_chunksize);
            while (pi[x]!=value) {
                if (jump and m_piv[x]==1) {
                    x  = s[piv_r1(x)]+(chunk*m_chunksize);
                    jump = false;
                } else {
                    x = pi[x]+(chunk*m_chunksize);
                }
            }
            return x_sls1(x+1)-x-(chunk*m_sigma)-1;
        }

        pair<size_type, value_type>	inverse_select(size_type i)const {
            uint64_t chunk = i/m_chunksize;
            bool jump=true;

            uint64_t x=i, value=i-(chunk*m_chunksize);
            while (pi[x]!=value) {
                if (jump and m_piv[x]==1) {
                    x  = s[piv_r1(x)]+(chunk*m_chunksize);
                    jump = false;
                } else {
                    x = pi[x]+(chunk*m_chunksize);
                }
            }
            uint64_t tmp = x_sls1(x+1);
            uint64_t c = tmp-x-(chunk*m_sigma)-1;
            uint64_t c_before_chunk;
            if (c==0) {
                if (chunk==0) {
                    c_before_chunk =0;
                } else {
                    c_before_chunk = b_sls0(chunk)-chunk+1;
                }
            } else {
                uint64_t ones_before_c = b_sls0(c*m_chunks)-(c*m_chunks)+1;
                c_before_chunk = b_sls0(c*m_chunks+chunk)-(c*m_chunks+chunk)+1-ones_before_c;
            }
            uint64_t c_in_chunk = tmp-x_sls0(c+1+chunk*m_sigma);
            return make_pair(c_before_chunk+c_in_chunk,c);
        }

        size_type select(size_type i, value_type c)const {

            uint64_t chunk, c_ones_before_chunk, pi_pos;

            if (c==0) {
                chunk = b_sls1(i)-i+1;
                if (chunk==0) {
                    c_ones_before_chunk = 0;
                } else {
                    c_ones_before_chunk = b_sls0(chunk)-chunk+1;
                }
                pi_pos = i-c_ones_before_chunk-1+chunk*m_chunksize;
            } else {
                uint64_t ones_before_c = b_sls0(c*m_chunks)-(c*m_chunks)+1;
                chunk = b_sls1(ones_before_c+i)-ones_before_c-(c*m_chunks)-i+1;
                c_ones_before_chunk = b_sls0(c*m_chunks+chunk)-(c*m_chunks+chunk)+1-ones_before_c;
                pi_pos = x_sls0(chunk*m_sigma+c+1)+(i-c_ones_before_chunk)-chunk*m_sigma-c-1;
            }

            return pi[pi_pos]+chunk*m_chunksize;
        }

        size_type rank(size_type i, value_type c)const {

            if (c>m_sigma-1) return 0;
            if (i<=0) return 0;

            uint64_t chunk = (i-1)/m_chunksize;
            uint64_t c_ones_before_chunk;

            if (c==0) {
                if (chunk==0) {
                    c_ones_before_chunk=0;
                } else {
                    c_ones_before_chunk=b_sls0(chunk)-chunk+1;
                }
            } else {
                uint64_t ones_before_c = b_sls0(c*m_chunks)-(c*m_chunks)+1;
                c_ones_before_chunk = b_sls0(c*m_chunks+chunk)-(c*m_chunks+chunk)+1-ones_before_c;
            }
            uint64_t c_ones_in_chunk = 0;

            size_type search_begin = x_sls0(chunk*m_sigma+1+c)-(chunk*m_sigma+1+c)+1;
            size_type search_end = x_sls0(chunk*m_sigma+2+c)-(chunk*m_sigma+2+c)+1;

            size_type val = (i-1)%m_chunksize;
            if (search_end-search_begin<50) { // After a short test, this seemd to be a good threshold
                while (search_begin < search_end and pi[search_begin] <= val) {
                    ++search_begin;
                    ++c_ones_in_chunk;
                }
            } else {
                c_ones_in_chunk = lower_bound(pi.begin()+search_begin, pi.begin()+search_end, val+1)-pi.begin()-search_begin;
            }
            return c_ones_before_chunk+c_ones_in_chunk;
        }

        //add members
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += write_member(m_chunks, out, child, "chunks");
            written_bytes += write_member(m_chunksize, out, child, "chunksize");
            written_bytes += write_member(m_t, out, child, "t");
            written_bytes += m_bv.serialize(out, child, "b");
            written_bytes += m_xv.serialize(out, child, "x");
            written_bytes += m_piv.serialize(out, child, "piv");
            written_bytes += pi.serialize(out, child, "pi");
            written_bytes += s.serialize(out, child, "s");
            written_bytes += b_sls0.serialize(out, child, "b_sls0");
            written_bytes += b_sls1.serialize(out, child, "b_sls1");
            written_bytes += x_sls1.serialize(out, child, "x_sls0");
            written_bytes += x_sls0.serialize(out, child, "x_sls1");
            written_bytes += piv_r1.serialize(out, child, "piv_r1");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};
