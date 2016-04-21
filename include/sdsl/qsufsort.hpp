/* qsufsort.c
   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.*/
/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file qsufsort.hpp
    \brief qsufsort.hpp contains the interface for the suffix array construction algorithm of Larsson.
	Larssons code was downloaded from http://www.larsson.dogma.net/qsufsort.c and adapted to the
	use of sdsl bitvectors.
	\author Simon Gog
*/

#ifndef INCLUDED_SDSL_QSUFSORT
#define INCLUDED_SDSL_QSUFSORT

#define DBG_OUT if(0)std::cout

#include <sdsl/int_vector.hpp>

namespace sdsl
{
namespace qsufsort
{

template<class int_vector_type=int_vector<> >
class sorter;

//void sort(int_iter text, int_iter sa, int64_t n, int64_t k, int64_t l);

//! Construct a suffix array for the sequence stored in a file.
/*!
 * \param sa		A reference to the resulting suffix array.
 * \param file  	Name of the file.
 * \param num_bytes Bytes per symbol in the file. I.e.
 *                  - num_bytes=1: byte sequence
 *                  - num_bytes=2: sequence of two byte symbols
 *                  - num_bytes=4: sequence of four byte symbols
 *                  - num_bytes=8: sequence of eight byte symbols.
 *                  - num_bytes=0: the algorithm assumes a serialized
 *                                 int_vector<> in the file and
 *                                 loads it.
 * \par Note
 *      If `int_vector_type` is `int_vector<>` then the bit-width of `sa` is
 *      the maximum of `bits::hi( max(sa.size()-1, 0) )` and the
 *      bit-width of the text.
 */
// TODO: problem when int_width==64!!!
template<class int_vector_type>
void construct_sa(int_vector_type& sa, const char* file, uint8_t num_bytes)
{
    sorter<int_vector_type> s;
    s.sort(sa, file, num_bytes);
}

template<class int_vector_type, class t_vec>
void construct_sa(int_vector_type& sa, t_vec& text)
{
    sorter<int_vector_type> s;
    s.sort(sa, text);
}

template<class int_vector_type>
class sorter
{
        typedef int_vector_type tIV;
        typedef typename tIV::iterator int_iter;
        typedef typename tIV::size_type size_type;

    private:

        int_iter m_SA,    // group array, ultimately suffix array.
                 m_VV;    // inverse array, ultimately inverse of SA.
        uint64_t m_rr,    // number of symbols aggregated by transform.
                 m_hh;    // length of already-sorted prefixes.
        uint8_t  m_msb;   // most significant bit position starting from 0
        uint64_t m_msb_mask;// mask for 1ULL<<msb

        inline int64_t to_sign(uint64_t x)const {
            return x & m_msb_mask ? -((int64_t)(x&~m_msb_mask)) : x;
        }
        // return the absolute value of integer x
        inline int64_t mark_pos(uint64_t x)const {
            return (x&~m_msb_mask);
        }
        // mark the number x as negative
        inline int64_t mark_neg(uint64_t x)const {
            return x|m_msb_mask;
        }
        // check if x is not negative
        inline bool     not_neg(uint64_t x)const {
            return !(x>>m_msb);
        }
        // check if x is negative
        inline bool     is_neg(uint64_t x)const {
            return x&m_msb_mask;
        }
        // returns the key of iterator p at the current sorting depth
        inline uint64_t key(const int_iter& p)const {
            return m_VV[*p+m_hh];
        }
        // swap the value of two iterators
        inline void swap(int_iter& p, int_iter& q)const {
            uint64_t tmp = *p; *p=*q; *q=tmp;
        }
        // select the median out of 3 elements
        inline const int_iter& med3(const int_iter& a, const int_iter& b, const int_iter& c)const {
            return key(a)<key(b)? (key(b)<key(c)?b:(key(a)<key(c)?c:a))
                       : (key(b)>key(c)?b:(key(a)>key(c)?c:a));
        }

        /* Subroutine for select_sort_split and sort_split. Sets group numbers for a
           group whose lowest position in m_SA is pl and highest position is pm.*/
        void update_group(int_iter pl, int_iter pm) {
            int64_t g=pm-m_SA;          /* group number.*/
            m_VV[*pl]=g;                /* update group number of first position.*/
            if (pl==pm)
                *pl=mark_neg(1);         /* one element, sorted group.*/
            else
                do                       /* more than one element, unsorted group.*/
                    m_VV[*++pl]=g;        /* update group numbers.*/
                while (pl<pm);
        }

        /* Quadratic sorting method to use for small subarrays. To be able to update
           group numbers consistently, a variant of selection sorting is used.*/
        void select_sort_split(const int_iter& p, int64_t n) {
            int_iter pa, pb, pi, pn;
            uint64_t f, v;

            pa=p;                        /* pa is start of group being picked out.*/
            pn=p+n-1;                    /* pn is last position of subarray.*/
            while (pa<pn) {
                for (pi=pb=(pa+1), f=key(pa); pi<=pn; ++pi)
                    if ((v=key(pi))<f) {
                        f=v;                /* f is smallest key found.*/
                        swap(pi, pa);       /* place smallest element at beginning.*/
                        pb=pa+1;            /* pb is position for elements equal to f.*/
                    } else if (v==f) {     /* if equal to smallest key.*/
                        swap(pi, pb);       /* place next to other smallest elements.*/
                        ++pb;
                    }
                update_group(pa, pb-1);   /* update group values for new group.*/
                pa=pb;                    /* continue sorting rest of the subarray.*/
            }
            if (pa==pn) {                /* check if last part is single element.*/
                m_VV[*pa]=pa-m_SA;
                *pa=mark_neg(1);          /* sorted group.*/
            }
        }

        /* Subroutine for sort_split, algorithm by Bentley & McIlroy.*/
        uint64_t choose_pivot(const int_iter& p, int64_t n) {
            int_iter pl, pm, pn;
            int64_t s;

            pm=p+(n>>1);                 /* small arrays, middle element.*/
            if (n>7LL) {
                pl=p;
                pn=p+n-1;
                if (n>40LL) {               /* big arrays, pseudomedian of 9.*/
                    s=n>>3;
                    pl=med3(pl, pl+s, pl+s+s);
                    pm=med3(pm-s, pm, pm+s);
                    pn=med3(pn-s-s, pn-s, pn);
                }
                pm=med3(pl, pm, pn);      /* midsize arrays, median of 3.*/
            }
            return key(pm);
        }

        /* Sorting routine called for each unsorted group. Sorts the array of integers
           (suffix numbers) of length n starting at p. The algorithm is a ternary-split
           quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
           Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
           function is based on Program 7.*/
        void sort_split(const int_iter& p, int64_t n) {
            int_iter pa, pb, pc, pd, pl, pm, pn;
            uint64_t f, v;
            int64_t  s, t;

            if (n<7) {                   /* multi-selection sort smallest arrays.*/
                select_sort_split(p, n);
                return;
            }

            v=choose_pivot(p, n);
            pa=pb=p;        // pa: iterator for equal elements
            pc=pd=p+n-1;    // pc = right border (inclusive)
            while (1) {                  /* split-end partition.*/
                while (pb<=pc && (f=key(pb))<=v) { // go to the right as long as there are keys smaller equal than v
                    if (f==v) {
                        swap(pa, pb); ++pa; // swap equal chars to the left
                    }
                    ++pb;
                }
                while (pc>=pb && (f=key(pc))>=v) { // go to the left as long as there are keys bigger or equal to v
                    if (f==v) {
                        swap(pc, pd); --pd; // swap equal chars to the right end
                    }
                    --pc;
                }
                if (pb>pc)
                    break;
                swap(pb, pc); // swap element > v (pb) to the third part and element < v (pc) to the second
                ++pb;
                --pc;
            }
            pn=p+n;
            if ((s=pa-p)>(t=pb-pa))
                s=t;
            for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
                swap(pl, pm);
            if ((s=pd-pc)>(t=pn-pd-1))
                s=t;
            for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
                swap(pl, pm);
            s=pb-pa;
            t=pd-pc;
            if (pa > pb) {
                if (s>0) {
                    std::cout<<"s="<<s<<">0 but should be <0; n="<<n<<std::endl;
                }
            }
            if (pc > pd) {
                if (t>0) {
                    std::cout<<"t="<<t<<">0 but should be <0; n="<<n<<std::endl;
                }
            }
            if (s>0)
                sort_split(p, s);
            update_group(p+s, p+n-t-1);
            if (t>0)
                sort_split(p+n-t, t);
        }

        /* Bucketsort for first iteration.

           Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
           at least once. x[n] is 0. (This is the corresponding output of transform.) k
           must be at most n+1. p is array of size n+1 whose contents are disregarded.

           Output: x is m_VV and p is m_SA after the initial sorting stage of the refined
           suffix sorting algorithm.*/
        void bucketsort(const int_iter& x, const int_iter& p, int64_t n, int64_t k) {
            int_iter pi;
            int64_t i, d, g;
            uint64_t c;

            for (pi=p; pi<p+k; ++pi)
                *pi=mark_neg(1);          /* mark linked lists empty.*/
            for (i=0; i<=n; ++i) {
                x[i]=p[c=x[i]];           /* insert in linked list.*/
                p[c]=i;
            }
            for (pi=p+k-1, i=n; pi>=p; --pi) {
                d=x[c=*pi];               /* c is position, d is next in list.*/
                x[c]=g=i;                 /* last position equals group number.*/
                if (not_neg(d)) {         /* if more than one element in group.*/
                    p[i--]=c;              /* p is permutation for the sorted x.*/
                    do {
                        d=x[c=d];           /* next in linked list.*/
                        x[c]=g;             /* group number in x.*/
                        p[i--]=c;           /* permutation in p.*/
                    } while (not_neg(d));
                } else
                    p[i--]=mark_neg(1);    /* one element, sorted group.*/
            }
        }

    public:

        /* Transforms the alphabet of x by attempting to aggregate several symbols into
           one, while preserving the suffix order of x. The alphabet may also be
           compacted, so that x on output comprises all integers of the new alphabet
           with no skipped numbers.

           Input: x is an array of size n+1 whose first n elements are positive
           integers in the range l...k-1. p is array of size n+1, used for temporary
           storage. q controls aggregation and compaction by defining the maximum value
           for any symbol during transformation: q must be at least k-l; if q<=n,
           compaction is guaranteed; if k-l>n, compaction is never done; if q is
           INT_MAX, the maximum number of symbols are aggregated into one.

           Output: Returns an integer j in the range 1...q representing the size of the
           new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
           set to the number of old symbols grouped into one. Only x[n] is 0.*/

        int64_t transform(const int_iter& x, const int_iter& p, int64_t n, int64_t k, int64_t l, int64_t q) {
            if (!(q >= k-l)) {
                std::cout << "q="<<q<<" k-l="<<k-l<<std::endl;
            }
            assert(q >= k-l);
            DBG_OUT<<"transform(n="<<n<<", k="<<k<<", l="<<l<<", q="<<q<<")"<<std::endl;
            uint64_t bb, cc,dd;
            int64_t jj;
            int_iter pi, pj;
            int s = bits::hi(k-l)+(k>l); /* s is number of bits in old symbol.*/
            uint8_t len = 0;                /* len is for overflow checking.*/
            m_rr = 0;
            for (bb=dd=0; (int)m_rr<n && (int)len < m_msb+1-s && (int64_t)(cc=dd<<s|(k-l)) <= q; ++m_rr, len+=s) {
                bb=bb<<s|(x[m_rr]-l+1);        /* bb is start of x in chunk alphabet.*/
                dd=cc;                      /* dd is max symbol in chunk alphabet.*/
            }
            DBG_OUT<<"m_rr="<<m_rr<<std::endl;
            uint64_t mm=(1ULL<<(m_rr-1)*s)-1;            /* mm masks off top old symbol from chunk.*/
            x[n]=l-1;                    /* emulate zero terminator.*/
            if ((int64_t)dd <= n) {                  /* if bucketing possible, compact alphabet.*/
                for (pi=p; pi<=p+dd; ++pi)
                    *pi=0;                 /* zero transformation table.*/
                for (pi=x+m_rr, cc=bb; pi<=x+n; ++pi) {
                    p[cc]=1;                /* mark used chunk symbol.*/
                    cc=(cc&mm)<<s|(*pi-l+1);  /* shift in next old symbol in chunk.*/
                }
                for (uint64_t i=1; i<m_rr; ++i) {     /* handle last r-1 positions.*/
                    p[cc]=1;                /* mark used chunk symbol.*/
                    cc=(cc&mm)<<s;            /* shift in next old symbol in chunk.*/
                }
                for (pi=p, jj=1; pi<=p+dd; ++pi)
                    if (*pi)
                        *pi=jj++;            /* j is new alphabet size.*/
                for (pi=x, pj=x+m_rr, cc=bb; pj<=x+n; ++pi, ++pj) {
                    *pi=p[cc];              /* transform to new alphabet.*/
                    cc=(cc&mm)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
                }
                while (pi<x+n) {          /* handle last r-1 positions.*/
                    *pi++=p[cc];            /* transform to new alphabet.*/
                    cc=(cc&mm)<<s;            /* shift right-end zero in chunk.*/
                }
            } else {                     /* bucketing not possible, don't compact.*/
                for (pi=x, pj=x+m_rr, cc=bb; pj<=x+n; ++pi, ++pj) {
                    *pi=cc;                 /* transform to new alphabet.*/
                    cc=(cc&mm)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
                }
                while (pi<x+n) {          /* handle last r-1 positions.*/
                    *pi++=cc;               /* transform to new alphabet.*/
                    cc=(cc&mm)<<s;            /* shift right-end zero in chunk.*/
                }
                jj=dd+1;                    /* new alphabet size.*/
            }
            x[n]=0;                      /* end-of-string symbol is zero.*/
            DBG_OUT<<"end transformation jj="<<jj<<std::endl;
            return jj;                    /* return new alphabet size.*/
        }

        /* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
           n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
           contents of x[n] is disregarded, the n-th symbol being regarded as
           end-of-string smaller than all other symbols.*/
        void sort(const int_iter& x, const int_iter& p, int64_t n, int64_t k, int64_t l) {
            int_iter pi, pk;
            m_VV=x;                         /* set global values.*/
            m_SA=p;
            if (n>=k-l) {                /* if bucketing possible,*/
                int64_t j = transform(m_VV, m_SA, n, k, l, n);
                DBG_OUT<<"begin bucketsort j="<<j<<std::endl;
                bucketsort(m_VV, m_SA, n, j);   /* bucketsort on first r positions.*/
                DBG_OUT<<"end bucketsort"<<std::endl;
            } else {
                transform(m_VV, m_SA, n, k, l, m_msb_mask-1);
                DBG_OUT<<"initialize SA begin"<<std::endl;
                for (int64_t i=0; i<=n; ++i)
                    m_SA[i]=i;                /* initialize I with suffix numbers.*/
                DBG_OUT<<"initialize SA end"<<std::endl;
                m_hh=0;
                sort_split(m_SA, n+1);       /* quicksort on first r positions.*/
            }
            m_hh=m_rr;                 /* number of symbols aggregated by transform.*/
//            while ( is_neg(*m_SA) and mark_pos(*m_SA) <= n) {
            while (to_sign(*m_SA) >= -n) {
//std::cout<<"m_hh="<<m_hh<<std::endl;
                DBG_OUT<<"SA = ";
//for(size_t iii=0; iii<=(size_t)n; ++iii){
//	uint64_t D = *(m_SA+iii);
//	printf("%c%lld ", is_neg(D)?'-':' ', mark_pos(D));
//}
                DBG_OUT<<std::endl;
                DBG_OUT<<"TEXT = ";
//for(size_t iii=0; iii<=(size_t)n; ++iii){
//	uint64_t D = *(m_VV+iii);
//	printf("%c%lld ", is_neg(D)?'-':' ', mark_pos(D));
//}
                DBG_OUT<<std::endl;
                DBG_OUT<<"*m_SA="<< to_sign(*m_SA) <<std::endl;
//std::cout<<"mark_pos(*m_SA)="<<mark_pos(*m_SA)<<std::endl;
                pi=m_SA;                     /* pi is first position of group.*/
                int64_t sl=0;              /* sl is length of sorted groups.*/
                DBG_OUT<<"m_hh="<<m_hh<<std::endl;
                do {
                    uint64_t s = *pi;
                    if (to_sign(s) < (int64_t)0) {
                        pi += mark_pos(s);   /* skip over sorted group.*/
                        sl += mark_pos(s);   /* add length to sl.*/
                    } else {
                        if (sl) {
                            *(pi-sl)=mark_neg(sl);     /* combine sorted groups before pi.*/
                            sl=0;
                        }
                        pk=m_SA+m_VV[s]+1;        /* pk-1 is last position of unsorted group.*/
                        sort_split(pi, pk-pi);
                        pi=pk;              /* next group.*/
                    }
                } while ((pi-m_SA) <= n);
                if (sl)                   /* if the array ends with a sorted group.*/
                    *(pi-sl)=mark_neg(sl);           /* combine sorted groups at end of m_SA.*/
                m_hh=2*m_hh;                    /* double sorted-depth.*/
                DBG_OUT<<"m_hh="<<m_hh<<std::endl;
            }
            for (int64_t i=0; i<=n; ++i) {        /* reconstruct suffix array from inverse.*/
                m_SA[m_VV[i]]=i;
            }
        }

        void do_sort(tIV& sa, tIV& x) {
            assert(x.size()>0);
            DBG_OUT<<"x.width()="<< (int)x.width() <<std::endl;
            DBG_OUT<<"x.size()="<<x.size()<<std::endl;
            DBG_OUT<<"sa.width()="<<(int)sa.width()<<std::endl;
            DBG_OUT<<"sa.size()="<<sa.size()<<std::endl;
            if (x.size() == 1) {
                sa = tIV(1, 0);
                return;
            }

            int64_t max_symbol = 0, min_symbol = x.width() < 64 ? bits::lo_set[x.width()] : 0x7FFFFFFFFFFFFFFFLL;

            for (size_type i=0; i < x.size()-1; ++i) {
                max_symbol = std::max(max_symbol, (int64_t)x[i]);
                min_symbol = std::min(min_symbol, (int64_t)x[i]);
            }

            if (0 == min_symbol) {
                throw std::logic_error("Text contains 0-symbol. Suffix array can not be constructed.");
            }
            if (x[x.size()-1] > 0) {
                throw std::logic_error("Last symbol is not 0-symbol. Suffix array can not be constructed.");
            }
            DBG_OUT<<"sorter: min_symbol="<<min_symbol<<std::endl;
            DBG_OUT<<"sorter: max_symbol="<<max_symbol<<std::endl;

            int64_t n = x.size()-1;
            DBG_OUT<<"x.size()-1="<<x.size()-1<<" n="<<n<<std::endl;
            uint8_t width = std::max(bits::hi(max_symbol)+2, bits::hi(n+1)+2);
            DBG_OUT<<"sorter: width="<<(int)width<<" max_symbol_width="<<bits::hi(max_symbol)+1<<" n_width="<< bits::hi(n) <<std::endl;
            util::expand_width(x, width);
            sa = x;
            if (sa.width() < x.width()) {
                throw std::logic_error("Fixed size suffix array is to small for the specified text.");
                return;
            }

            m_msb = sa.width()-1;
            m_msb_mask = 1ULL<<m_msb;
            DBG_OUT<<"sorter: m_msb="<< (int)m_msb <<" m_msb_mask="<<m_msb_mask<<std::endl;
            sort(x.begin(), sa.begin(), x.size()-1, max_symbol+1, min_symbol);
        }


        void sort(tIV& sa, const char* file_name, uint8_t num_bytes) {
            DBG_OUT<<"sorter: sort("<<file_name<<")"<<std::endl;
            DBG_OUT<<"sizeof(int_vector<>::difference_type)="<<sizeof(int_vector<>::difference_type)<<std::endl;
            util::clear(sa); // free space for sa
            tIV x;
            if (num_bytes == 0 and typeid(typename tIV::reference) == typeid(uint64_t)) {
                DBG_OUT<<"sorter: use int_vector<64>"<<std::endl;
                int_vector<> temp;
                load_vector_from_file(temp, file_name, num_bytes);
                x.resize(temp.size());
                for (size_type i=0; i<temp.size(); ++i) x[i] = temp[i];
            } else {
                load_vector_from_file(x, file_name, num_bytes);
                util::bit_compress(x);
            }
            do_sort(sa, x);
        }

        template<class t_vec>
        void sort(tIV& sa, t_vec& text) {
            tIV x;
            x.resize(text.size());
            for (size_type i=0; i<text.size(); ++i) x[i] = text[i];
            do_sort(sa, x);
        }
};

} // end namespace qsufsort

} // end namespace sdsl

#endif
