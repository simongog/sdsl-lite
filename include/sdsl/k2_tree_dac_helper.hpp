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
/*! \file k2_treap_helper.hpp
    \brief k2_treap_helper.hpp contains helper functions and definitions for a k^2-treap implementation.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREE_BITS_HELPER
#define INCLUDED_SDSL_K2_TREE_BITS_HELPER

#include "vectors.hpp"
#include "bits.hpp"
#include <tuple>
#include <algorithm>
#include <iterator>
#include <vector>
#include <complex>
#include <queue>
#include <array>


// Basics

// #include "basics.h" included later to avoid macro recursion for malloc
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>

//! Namespace for the succinct data structure library.
namespace sdsl {
    #define mask31 0x0000001F

    /*numero de bits del entero de la maquina*/
    #define W 32
    /* W-1 */
    #define Wminusone 31
    /*numero de bits del entero de la maquina*/
    #define WW 64
    /*bits para hacer la mascara para contar mas rapido*/
    #define bitsM 8
    /*bytes que hacen una palabra */
    #define BW 4
    #ifndef uint
    #define uint unsigned int
    #endif
    #ifndef ulong
    #define ulong unsigned long
    #endif
    #define size_uchar 256

    #ifndef byte
    #define byte unsigned char
    #endif

    /* reads bit p from e */
    #define bitgetchar(e, p) ((((e)[(p)/bitsM] >> ((p)%bitsM))) & 1)
    /* sets bit p in e */
    #define bitsetchar(e, p) ((e)[(p)/bitsM] |= (1<<((p)%bitsM)))
    /* cleans bit p in e */
    #define bitcleanchar(e, p) ((e)[(p)/bitsM] &= ~(1<<((p)%bitsM)))


    /* reads bit p from e */
    #define bitget(e, p) ((((e)[(p)/W] >> ((p)%W))) & 1)
    /* sets bit p in e */
    #define bitset(e, p) ((e)[(p)/W] |= (1<<((p)%W)))
    /* cleans bit p in e */
    #define bitclean(e, p) ((e)[(p)/W] &= ~(1<<((p)%W)))


    typedef int_vector<>::size_type size_type;
    typedef unsigned char uchar;
/** Number of bits in a byte */
    const uint kByteBits = 8;
/** Number of bits in an unsigned char */
    const uint kUcharBits = kByteBits * sizeof(uchar);

    namespace k2_dac_helper {

        const unsigned char __popcount_tab[] = {
                0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
        };

        const unsigned char select_tab[] = {
                0, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                6, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                7, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                6, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                8, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                6, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                7, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
                6, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1,
        };

        const unsigned char prev_tab[] = {
                0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        };

        // Memory management

        void *Malloc(int n) {
            void *p;
            if (n == 0) return NULL;
            p = (void *) malloc(n);
            if (p == NULL) {
                fprintf(stderr, "Could not allocate %i bytes\n", n);
                exit(1);
            }
            return p;
        }

        void Free(void *p) {
            if (p) free(p);
        }

        void *Realloc(void *p, int n) {
            if (p == NULL) return Malloc(n);
            if (n == 0) {
                Free(p);
                return NULL;
            }
            p = (void *) realloc(p, n);
            if (p == NULL) {
                fprintf(stderr, "Could not allocate %i bytes\n", n);
                exit(1);
            }
            return p;
        }

        // bits needed to represent a number between 0 and n

        uint bits(uint n) {
            uint b = 0;
            while (n) {
                b++;
                n >>= 1;
            }
            return b;
        }

        // returns e[p..p+len-1], assuming len <= W

        uint bitread(uint *e, uint p, uint len) {
            uint answ;
            e += p / W;
            p %= W;
            answ = *e >> p;
            if (len == W) {
                if (p) answ |= (*(e + 1)) << (W - p);
            }
            else {
                if (p + len > W) answ |= (*(e + 1)) << (W - p);
                answ &= (1 << len) - 1;
            }
            return answ;
        }


        // writes e[p..p+len-1] = s, len <= W

        void bitwrite(register uint *e, register uint p,
                      register uint len, register uint s) {
            e += p / W;
            p %= W;
            if (len == W) {
                *e |= (*e & ((1 << p) - 1)) | (s << p);
                if (!p) return;
                e++;
                *e = (*e & ~((1 << p) - 1)) | (s >> (W - p));
            }
            else {
                if (p + len <= W) {
                    *e = (*e & ~(((1u << len) - 1) << p)) | (s << p);
                    return;
                }
                *e = (*e & ((1 << p) - 1)) | (s << p);
                e++;
                len -= W - p;
                *e = (*e & ~((1 << len) - 1)) | (s >> (W - p));
            }
        }
        // writes e[p..p+len-1] = 0

        void bitzero(register uint *e, register uint p,
                     register uint len) {
            e += p / W;
            p %= W;
            if (p + len >= W) {
                *e &= ~((1 << p) - 1);
                len -= p;
                e++;
                p = 0;
            }
            while (len >= W) {
                *e++ = 0;
                len -= W;
            }
            if (len > 0)
                *e &= ~(((1 << len) - 1) << p);
        }


        uint GetField(uint *A, register uint len, register uint index) {
            register uint i = index * len / W, j = index * len - W * i, result;
            if (j + len <= W)
                result = (A[i] << (W - j - len)) >> (W - len);
            else {
                result = A[i] >> j;
                result = result | (A[i + 1] << (WW - j - len)) >> (W - len);
            }
            return result;
        }


        void SetField(uint *A, register uint len, register uint index, register uint x) {
            uint i = index * len / W, j = index * len - i * W;
            uint mask = ((j + len) < W ? ~0u << (j + len) : 0) | ((W - j) < W ? ~0u >> (W - j) : 0);
            A[i] = (A[i] & mask) | x << j;
            if (j + len > W) {
                mask = ((~0u) << (len + j - W));
                A[i + 1] = (A[i + 1] & mask) | x >> (W - j);
            }
        }


        uint GetVarField(uint *A, register uint ini, register uint fin) {
            register uint i = ini / W, j = ini - W * i, result;
            register uint len = (fin - ini + 1);
            if (j + len <= W)
                result = (A[i] << (W - j - len)) >> (W - len);
            else {
                result = A[i] >> j;
                result = result | (A[i + 1] << (WW - j - len)) >> (W - len);
            }
            return result;
        }


        void SetVarField(uint *A, register uint ini, register uint fin, register uint x) {
            uint i = ini / W, j = ini - i * W;
            uint len = (fin - ini + 1);
            uint mask = ((j + len) < W ? ~0u << (j + len) : 0) | ((W - j) < W ? ~0u >> (W - j) : 0);
            A[i] = (A[i] & mask) | x << j;
            if (j + len > W) {
                mask = ((~0u) << (len + j - W));
                A[i + 1] = (A[i + 1] & mask) | x >> (W - j);
            }
        }


        unsigned GetFieldW32(uint *A, register uint index) {
            return A[index];
        }


        void SetField32(uint *A, register uint index, register uint x) {
            A[index] = x;
        }


        unsigned GetFieldW16(uint *A, register uint index) {
            register uint i = index / 2, j = (index & 1) << 4, result;
            result = (A[i] << (16 - j)) >> (16);
            return result;
        }

        unsigned GetFieldW4(uint *A, register uint index) {
            register uint i = index / 8, j = (index & 0x7) << 2;
            /*register unsigned i=index/8, j=index*4-32*i; */
            return (A[i] << (28 - j)) >> (28);
        }


        uint popcount(register int x) {
            return __popcount_tab[(x >> 0) & 0xff] + __popcount_tab[(x >> 8) & 0xff] +
                   __popcount_tab[(x >> 16) & 0xff] + __popcount_tab[(x >> 24) & 0xff];
        }


        uint popcount16(register int x) {
            return __popcount_tab[x & 0xff] + __popcount_tab[(x >> 8) & 0xff];
        }


        uint popcount8(register int x) {
            return __popcount_tab[x & 0xff];
        }


    }

}

#endif