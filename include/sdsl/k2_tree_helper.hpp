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
#ifndef INCLUDED_SDSL_K2_TREE_HELPER
#define INCLUDED_SDSL_K2_TREE_HELPER

#include "vectors.hpp"
#include "bits.hpp"
#include <tuple>
#include <algorithm>
#include <iterator>
#include <vector>
#include <complex>
#include <queue>
#include <array>
#include <bitset>

//! Namespace for the succinct data structure library.
namespace sdsl {
    typedef int_vector<>::size_type size_type;

    namespace k2_treap_ns {

        template<typename t_x=uint64_t, typename t_y=uint64_t>
        std::vector<std::pair<t_x, t_y>> read(std::vector<int_vector_buffer<> * > &bufs) {
            typedef std::vector<std::pair<t_x, t_y>> t_tuple_vec;
            t_tuple_vec v = t_tuple_vec(bufs[0]->size());
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<0>(v[j]) = (*(bufs[0]))[j];
            }
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<1>(v[j]) = (*(bufs[1]))[j];
            }

            return v;
        }

        /**
        * Calculates the smalles integer gretear or equal to x/y
        */
        template<typename T>
        inline T div_ceil(T x, T y) {
            static_assert(std::is_integral<T>::value, "Parameter is not integral type");
            return (x % y) ? x / y + 1 : x / y;
        }

        //FIXME: assumes one Byte has 8 Bits, which is not always the case according to c++ standard
        static const uint kUcharBits = 8 * sizeof(unsigned char);

        // Precomputed value for fast k^2 treap operations
        template<uint8_t t_k>
        struct precomp {
            static struct impl {
                uint64_t exp[65];

                impl() {
                    exp[0] = 1;
                    for (uint8_t i = 1; i < 65; ++i) {
                        exp[i] = t_k * exp[i - 1];
                    }
                }
            } data;

            static uint64_t exp(uint8_t l) {
                return data.exp[l];
            }

            static uint64_t divexp(uint64_t x, uint8_t l) {
                return x / data.exp[l];
            }

            static uint64_t modexp(uint64_t x, uint8_t l) {
                return x % data.exp[l];
            }
        };

        template<>
        struct precomp<2> {
            static uint64_t exp(uint8_t l) {
                return 1ULL << l;
            }

            static uint64_t divexp(uint64_t x, uint8_t l) {
                return x >> l;
            }

            static uint64_t modexp(uint64_t x, uint8_t l) {
                return x & bits::lo_set[l];
            }
        };

        template<>
        struct precomp<4> {
            static uint64_t exp(uint8_t l) {
                return 1ULL << (2 * l);
            }

            static uint64_t divexp(uint64_t x, uint8_t l) {
                return x >> (2 * l);
            }

            static uint64_t modexp(uint64_t x, uint8_t l) {
                return x & bits::lo_set[2 * l];
            }
        };

        template<>
        struct precomp<8> {
            static uint64_t exp(uint8_t l) {
                return 1ULL << (3 * l);
            }

            static uint64_t divexp(uint64_t x, uint8_t l) {
                return x >> (3 * l);
            }

            static uint64_t modexp(uint64_t x, uint8_t l) {
                return x & bits::lo_set[3 * l];
            }
        };

        template<>
        struct precomp<16> {
            static uint64_t exp(uint8_t l) {
                return 1ULL << (4 * l);
            }

            static uint64_t divexp(uint64_t x, uint8_t l) {
                return x >> (4 * l);
            }

            static uint64_t modexp(uint64_t x, uint8_t l) {
                return x & bits::lo_set[4 * l];
            }
        };


        template<uint8_t t_k>
        typename precomp<t_k>::impl precomp<t_k>::data;

        typedef std::complex<uint64_t> t_p;
        typedef t_p point_type;

        struct node_type {
            uint8_t t;   // level; size of node 1<<t
            t_p p;   // lower left corner; upper left in matrix!
            uint64_t idx; // index in bp

            node_type() = default;

            node_type(uint8_t _t, t_p _p, uint64_t _idx) : t(_t), p(_p), idx(_idx) {}

            node_type(node_type &&) = default;

            node_type(const node_type &) = default;

            node_type &operator=(node_type &&) = default;

            node_type &operator=(const node_type &) = default;
        };

        //TODO: add more preshifted Morton Tables

        static constexpr uint_fast16_t morton_x_256[256] =
                {
                        0, 1, 4, 5, 16, 17, 20, 21,
                        64, 65, 68, 69, 80, 81, 84, 85,
                        256, 257, 260, 261, 272, 273, 276, 277,
                        320, 321, 324, 325, 336, 337, 340, 341,
                        1024, 1025, 1028, 1029, 1040, 1041, 1044, 1045,
                        1088, 1089, 1092, 1093, 1104, 1105, 1108, 1109,
                        1280, 1281, 1284, 1285, 1296, 1297, 1300, 1301,
                        1344, 1345, 1348, 1349, 1360, 1361, 1364, 1365,
                        4096, 4097, 4100, 4101, 4112, 4113, 4116, 4117,
                        4160, 4161, 4164, 4165, 4176, 4177, 4180, 4181,
                        4352, 4353, 4356, 4357, 4368, 4369, 4372, 4373,
                        4416, 4417, 4420, 4421, 4432, 4433, 4436, 4437,
                        5120, 5121, 5124, 5125, 5136, 5137, 5140, 5141,
                        5184, 5185, 5188, 5189, 5200, 5201, 5204, 5205,
                        5376, 5377, 5380, 5381, 5392, 5393, 5396, 5397,
                        5440, 5441, 5444, 5445, 5456, 5457, 5460, 5461,
                        16384, 16385, 16388, 16389, 16400, 16401, 16404, 16405,
                        16448, 16449, 16452, 16453, 16464, 16465, 16468, 16469,
                        16640, 16641, 16644, 16645, 16656, 16657, 16660, 16661,
                        16704, 16705, 16708, 16709, 16720, 16721, 16724, 16725,
                        17408, 17409, 17412, 17413, 17424, 17425, 17428, 17429,
                        17472, 17473, 17476, 17477, 17488, 17489, 17492, 17493,
                        17664, 17665, 17668, 17669, 17680, 17681, 17684, 17685,
                        17728, 17729, 17732, 17733, 17744, 17745, 17748, 17749,
                        20480, 20481, 20484, 20485, 20496, 20497, 20500, 20501,
                        20544, 20545, 20548, 20549, 20560, 20561, 20564, 20565,
                        20736, 20737, 20740, 20741, 20752, 20753, 20756, 20757,
                        20800, 20801, 20804, 20805, 20816, 20817, 20820, 20821,
                        21504, 21505, 21508, 21509, 21520, 21521, 21524, 21525,
                        21568, 21569, 21572, 21573, 21584, 21585, 21588, 21589,
                        21760, 21761, 21764, 21765, 21776, 21777, 21780, 21781,
                        21824, 21825, 21828, 21829, 21840, 21841, 21844, 21845
                };

        static constexpr uint_fast16_t morton_y_256[256] =
                {
                        0, 2, 8, 10, 32, 34, 40, 42,
                        128, 130, 136, 138, 160, 162, 168, 170,
                        512, 514, 520, 522, 544, 546, 552, 554,
                        640, 642, 648, 650, 672, 674, 680, 682,
                        2048, 2050, 2056, 2058, 2080, 2082, 2088, 2090,
                        2176, 2178, 2184, 2186, 2208, 2210, 2216, 2218,
                        2560, 2562, 2568, 2570, 2592, 2594, 2600, 2602,
                        2688, 2690, 2696, 2698, 2720, 2722, 2728, 2730,
                        8192, 8194, 8200, 8202, 8224, 8226, 8232, 8234,
                        8320, 8322, 8328, 8330, 8352, 8354, 8360, 8362,
                        8704, 8706, 8712, 8714, 8736, 8738, 8744, 8746,
                        8832, 8834, 8840, 8842, 8864, 8866, 8872, 8874,
                        10240, 10242, 10248, 10250, 10272, 10274, 10280, 10282,
                        10368, 10370, 10376, 10378, 10400, 10402, 10408, 10410,
                        10752, 10754, 10760, 10762, 10784, 10786, 10792, 10794,
                        10880, 10882, 10888, 10890, 10912, 10914, 10920, 10922,
                        32768, 32770, 32776, 32778, 32800, 32802, 32808, 32810,
                        32896, 32898, 32904, 32906, 32928, 32930, 32936, 32938,
                        33280, 33282, 33288, 33290, 33312, 33314, 33320, 33322,
                        33408, 33410, 33416, 33418, 33440, 33442, 33448, 33450,
                        34816, 34818, 34824, 34826, 34848, 34850, 34856, 34858,
                        34944, 34946, 34952, 34954, 34976, 34978, 34984, 34986,
                        35328, 35330, 35336, 35338, 35360, 35362, 35368, 35370,
                        35456, 35458, 35464, 35466, 35488, 35490, 35496, 35498,
                        40960, 40962, 40968, 40970, 40992, 40994, 41000, 41002,
                        41088, 41090, 41096, 41098, 41120, 41122, 41128, 41130,
                        41472, 41474, 41480, 41482, 41504, 41506, 41512, 41514,
                        41600, 41602, 41608, 41610, 41632, 41634, 41640, 41642,
                        43008, 43010, 43016, 43018, 43040, 43042, 43048, 43050,
                        43136, 43138, 43144, 43146, 43168, 43170, 43176, 43178,
                        43520, 43522, 43528, 43530, 43552, 43554, 43560, 43562,
                        43648, 43650, 43656, 43658, 43680, 43682, 43688, 43690
                };


        static constexpr uint_fast16_t morton_2_bitwise_y_256[256] =
                {
                        0x0000, 0x0004, 0x0008, 0x000c, 0x0040, 0x0044, 0x0048, 0x004c, 0x0080, 0x0084, 0x0088, 0x008c,
                        0x00c0, 0x00c4, 0x00c8, 0x00cc,
                        0x0400, 0x0404, 0x0408, 0x040c, 0x0440, 0x0444, 0x0448, 0x044c, 0x0480, 0x0484, 0x0488, 0x048c,
                        0x04c0, 0x04c4, 0x04c8, 0x04cc,
                        0x0800, 0x0804, 0x0808, 0x080c, 0x0840, 0x0844, 0x0848, 0x084c, 0x0880, 0x0884, 0x0888, 0x088c,
                        0x08c0, 0x08c4, 0x08c8, 0x08cc,
                        0x0c00, 0x0c04, 0x0c08, 0x0c0c, 0x0c40, 0x0c44, 0x0c48, 0x0c4c, 0x0c80, 0x0c84, 0x0c88, 0x0c8c,
                        0x0cc0, 0x0cc4, 0x0cc8, 0x0ccc,
                        0x4000, 0x4004, 0x4008, 0x400c, 0x4040, 0x4044, 0x4048, 0x404c, 0x4080, 0x4084, 0x4088, 0x408c,
                        0x40c0, 0x40c4, 0x40c8, 0x40cc,
                        0x4400, 0x4404, 0x4408, 0x440c, 0x4440, 0x4444, 0x4448, 0x444c, 0x4480, 0x4484, 0x4488, 0x448c,
                        0x44c0, 0x44c4, 0x44c8, 0x44cc,
                        0x4800, 0x4804, 0x4808, 0x480c, 0x4840, 0x4844, 0x4848, 0x484c, 0x4880, 0x4884, 0x4888, 0x488c,
                        0x48c0, 0x48c4, 0x48c8, 0x48cc,
                        0x4c00, 0x4c04, 0x4c08, 0x4c0c, 0x4c40, 0x4c44, 0x4c48, 0x4c4c, 0x4c80, 0x4c84, 0x4c88, 0x4c8c,
                        0x4cc0, 0x4cc4, 0x4cc8, 0x4ccc,
                        0x8000, 0x8004, 0x8008, 0x800c, 0x8040, 0x8044, 0x8048, 0x804c, 0x8080, 0x8084, 0x8088, 0x808c,
                        0x80c0, 0x80c4, 0x80c8, 0x80cc,
                        0x8400, 0x8404, 0x8408, 0x840c, 0x8440, 0x8444, 0x8448, 0x844c, 0x8480, 0x8484, 0x8488, 0x848c,
                        0x84c0, 0x84c4, 0x84c8, 0x84cc,
                        0x8800, 0x8804, 0x8808, 0x880c, 0x8840, 0x8844, 0x8848, 0x884c, 0x8880, 0x8884, 0x8888, 0x888c,
                        0x88c0, 0x88c4, 0x88c8, 0x88cc,
                        0x8c00, 0x8c04, 0x8c08, 0x8c0c, 0x8c40, 0x8c44, 0x8c48, 0x8c4c, 0x8c80, 0x8c84, 0x8c88, 0x8c8c,
                        0x8cc0, 0x8cc4, 0x8cc8, 0x8ccc,
                        0xc000, 0xc004, 0xc008, 0xc00c, 0xc040, 0xc044, 0xc048, 0xc04c, 0xc080, 0xc084, 0xc088, 0xc08c,
                        0xc0c0, 0xc0c4, 0xc0c8, 0xc0cc,
                        0xc400, 0xc404, 0xc408, 0xc40c, 0xc440, 0xc444, 0xc448, 0xc44c, 0xc480, 0xc484, 0xc488, 0xc48c,
                        0xc4c0, 0xc4c4, 0xc4c8, 0xc4cc,
                        0xc800, 0xc804, 0xc808, 0xc80c, 0xc840, 0xc844, 0xc848, 0xc84c, 0xc880, 0xc884, 0xc888, 0xc88c,
                        0xc8c0, 0xc8c4, 0xc8c8, 0xc8cc,
                        0xcc00, 0xcc04, 0xcc08, 0xcc0c, 0xcc40, 0xcc44, 0xcc48, 0xcc4c, 0xcc80, 0xcc84, 0xcc88, 0xcc8c,
                        0xccc0, 0xccc4, 0xccc8, 0xcccc
                };


        static constexpr uint_fast16_t morton_2_bitwise_x_256[256] = {
                0x0, 0x1, 0x2, 0x3, 0x10, 0x11, 0x12, 0x13, 0x20, 0x21, 0x22, 0x23, 0x30, 0x31, 0x32, 0x33, 0x100,
                0x101, 0x102, 0x103, 0x110, 0x111, 0x112, 0x113, 0x120, 0x121, 0x122, 0x123, 0x130, 0x131, 0x132, 0x133,
                0x200, 0x201, 0x202, 0x203, 0x210, 0x211, 0x212, 0x213, 0x220, 0x221, 0x222, 0x223, 0x230, 0x231, 0x232,
                0x233, 0x300, 0x301, 0x302, 0x303, 0x310, 0x311, 0x312, 0x313, 0x320, 0x321, 0x322, 0x323, 0x330, 0x331,
                0x332, 0x333, 0x1000, 0x1001, 0x1002, 0x1003, 0x1010, 0x1011, 0x1012, 0x1013, 0x1020, 0x1021, 0x1022,
                0x1023, 0x1030, 0x1031, 0x1032, 0x1033, 0x1100, 0x1101, 0x1102, 0x1103, 0x1110, 0x1111, 0x1112, 0x1113,
                0x1120, 0x1121, 0x1122, 0x1123, 0x1130, 0x1131, 0x1132, 0x1133, 0x1200, 0x1201, 0x1202, 0x1203, 0x1210,
                0x1211, 0x1212, 0x1213, 0x1220, 0x1221, 0x1222, 0x1223, 0x1230, 0x1231, 0x1232, 0x1233, 0x1300, 0x1301,
                0x1302, 0x1303, 0x1310, 0x1311, 0x1312, 0x1313, 0x1320, 0x1321, 0x1322, 0x1323, 0x1330, 0x1331, 0x1332,
                0x1333, 0x2000, 0x2001, 0x2002, 0x2003, 0x2010, 0x2011, 0x2012, 0x2013, 0x2020, 0x2021, 0x2022, 0x2023,
                0x2030, 0x2031, 0x2032, 0x2033, 0x2100, 0x2101, 0x2102, 0x2103, 0x2110, 0x2111, 0x2112, 0x2113, 0x2120,
                0x2121, 0x2122, 0x2123, 0x2130, 0x2131, 0x2132, 0x2133, 0x2200, 0x2201, 0x2202, 0x2203, 0x2210, 0x2211,
                0x2212, 0x2213, 0x2220, 0x2221, 0x2222, 0x2223, 0x2230, 0x2231, 0x2232, 0x2233, 0x2300, 0x2301, 0x2302,
                0x2303, 0x2310, 0x2311, 0x2312, 0x2313, 0x2320, 0x2321, 0x2322, 0x2323, 0x2330, 0x2331, 0x2332, 0x2333,
                0x3000, 0x3001, 0x3002, 0x3003, 0x3010, 0x3011, 0x3012, 0x3013, 0x3020, 0x3021, 0x3022, 0x3023, 0x3030,
                0x3031, 0x3032, 0x3033, 0x3100, 0x3101, 0x3102, 0x3103, 0x3110, 0x3111, 0x3112, 0x3113, 0x3120, 0x3121,
                0x3122, 0x3123, 0x3130, 0x3131, 0x3132, 0x3133, 0x3200, 0x3201, 0x3202, 0x3203, 0x3210, 0x3211, 0x3212,
                0x3213, 0x3220, 0x3221, 0x3222, 0x3223, 0x3230, 0x3231, 0x3232, 0x3233, 0x3300, 0x3301, 0x3302, 0x3303,
                0x3310, 0x3311, 0x3312, 0x3313, 0x3320, 0x3321, 0x3322, 0x3323, 0x3330, 0x3331, 0x3332, 0x3333
        };

        static constexpr uint_fast16_t morton_3_bitwise_x_512[512] =
                {
                        0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x80,
                        0x81,
                        0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7, 0x100,
                        0x101,
                        0x102, 0x103, 0x104, 0x105, 0x106, 0x107, 0x140, 0x141, 0x142, 0x143, 0x144, 0x145, 0x146,
                        0x147,
                        0x180, 0x181, 0x182, 0x183, 0x184, 0x185, 0x186, 0x187, 0x1c0, 0x1c1, 0x1c2, 0x1c3, 0x1c4,
                        0x1c5,
                        0x1c6, 0x1c7, 0x1000, 0x1001, 0x1002, 0x1003, 0x1004, 0x1005, 0x1006, 0x1007, 0x1040, 0x1041,
                        0x1042, 0x1043, 0x1044, 0x1045, 0x1046, 0x1047, 0x1080, 0x1081, 0x1082, 0x1083, 0x1084, 0x1085,
                        0x1086, 0x1087, 0x10c0, 0x10c1, 0x10c2, 0x10c3, 0x10c4, 0x10c5, 0x10c6, 0x10c7, 0x1100, 0x1101,
                        0x1102, 0x1103, 0x1104, 0x1105, 0x1106, 0x1107, 0x1140, 0x1141, 0x1142, 0x1143, 0x1144, 0x1145,
                        0x1146, 0x1147, 0x1180, 0x1181, 0x1182, 0x1183, 0x1184, 0x1185, 0x1186, 0x1187, 0x11c0, 0x11c1,
                        0x11c2, 0x11c3, 0x11c4, 0x11c5, 0x11c6, 0x11c7, 0x2000, 0x2001, 0x2002, 0x2003, 0x2004, 0x2005,
                        0x2006, 0x2007, 0x2040, 0x2041, 0x2042, 0x2043, 0x2044, 0x2045, 0x2046, 0x2047, 0x2080, 0x2081,
                        0x2082, 0x2083, 0x2084, 0x2085, 0x2086, 0x2087, 0x20c0, 0x20c1, 0x20c2, 0x20c3, 0x20c4, 0x20c5,
                        0x20c6, 0x20c7, 0x2100, 0x2101, 0x2102, 0x2103, 0x2104, 0x2105, 0x2106, 0x2107, 0x2140, 0x2141,
                        0x2142, 0x2143, 0x2144, 0x2145, 0x2146, 0x2147, 0x2180, 0x2181, 0x2182, 0x2183, 0x2184, 0x2185,
                        0x2186, 0x2187, 0x21c0, 0x21c1, 0x21c2, 0x21c3, 0x21c4, 0x21c5, 0x21c6, 0x21c7, 0x3000, 0x3001,
                        0x3002, 0x3003, 0x3004, 0x3005, 0x3006, 0x3007, 0x3040, 0x3041, 0x3042, 0x3043, 0x3044, 0x3045,
                        0x3046, 0x3047, 0x3080, 0x3081, 0x3082, 0x3083, 0x3084, 0x3085, 0x3086, 0x3087, 0x30c0, 0x30c1,
                        0x30c2, 0x30c3, 0x30c4, 0x30c5, 0x30c6, 0x30c7, 0x3100, 0x3101, 0x3102, 0x3103, 0x3104, 0x3105,
                        0x3106, 0x3107, 0x3140, 0x3141, 0x3142, 0x3143, 0x3144, 0x3145, 0x3146, 0x3147, 0x3180, 0x3181,
                        0x3182, 0x3183, 0x3184, 0x3185, 0x3186, 0x3187, 0x31c0, 0x31c1, 0x31c2, 0x31c3, 0x31c4, 0x31c5,
                        0x31c6, 0x31c7, 0x4000, 0x4001, 0x4002, 0x4003, 0x4004, 0x4005, 0x4006, 0x4007, 0x4040, 0x4041,
                        0x4042, 0x4043, 0x4044, 0x4045, 0x4046, 0x4047, 0x4080, 0x4081, 0x4082, 0x4083, 0x4084, 0x4085,
                        0x4086, 0x4087, 0x40c0, 0x40c1, 0x40c2, 0x40c3, 0x40c4, 0x40c5, 0x40c6, 0x40c7, 0x4100, 0x4101,
                        0x4102, 0x4103, 0x4104, 0x4105, 0x4106, 0x4107, 0x4140, 0x4141, 0x4142, 0x4143, 0x4144, 0x4145,
                        0x4146, 0x4147, 0x4180, 0x4181, 0x4182, 0x4183, 0x4184, 0x4185, 0x4186, 0x4187, 0x41c0, 0x41c1,
                        0x41c2, 0x41c3, 0x41c4, 0x41c5, 0x41c6, 0x41c7, 0x5000, 0x5001, 0x5002, 0x5003, 0x5004, 0x5005,
                        0x5006, 0x5007, 0x5040, 0x5041, 0x5042, 0x5043, 0x5044, 0x5045, 0x5046, 0x5047, 0x5080, 0x5081,
                        0x5082, 0x5083, 0x5084, 0x5085, 0x5086, 0x5087, 0x50c0, 0x50c1, 0x50c2, 0x50c3, 0x50c4, 0x50c5,
                        0x50c6, 0x50c7, 0x5100, 0x5101, 0x5102, 0x5103, 0x5104, 0x5105, 0x5106, 0x5107, 0x5140, 0x5141,
                        0x5142, 0x5143, 0x5144, 0x5145, 0x5146, 0x5147, 0x5180, 0x5181, 0x5182, 0x5183, 0x5184, 0x5185,
                        0x5186, 0x5187, 0x51c0, 0x51c1, 0x51c2, 0x51c3, 0x51c4, 0x51c5, 0x51c6, 0x51c7, 0x6000, 0x6001,
                        0x6002, 0x6003, 0x6004, 0x6005, 0x6006, 0x6007, 0x6040, 0x6041, 0x6042, 0x6043, 0x6044, 0x6045,
                        0x6046, 0x6047, 0x6080, 0x6081, 0x6082, 0x6083, 0x6084, 0x6085, 0x6086, 0x6087, 0x60c0, 0x60c1,
                        0x60c2, 0x60c3, 0x60c4, 0x60c5, 0x60c6, 0x60c7, 0x6100, 0x6101, 0x6102, 0x6103, 0x6104, 0x6105,
                        0x6106, 0x6107, 0x6140, 0x6141, 0x6142, 0x6143, 0x6144, 0x6145, 0x6146, 0x6147, 0x6180, 0x6181,
                        0x6182, 0x6183, 0x6184, 0x6185, 0x6186, 0x6187, 0x61c0, 0x61c1, 0x61c2, 0x61c3, 0x61c4, 0x61c5,
                        0x61c6, 0x61c7, 0x7000, 0x7001, 0x7002, 0x7003, 0x7004, 0x7005, 0x7006, 0x7007, 0x7040, 0x7041,
                        0x7042, 0x7043, 0x7044, 0x7045, 0x7046, 0x7047, 0x7080, 0x7081, 0x7082, 0x7083, 0x7084, 0x7085,
                        0x7086, 0x7087, 0x70c0, 0x70c1, 0x70c2, 0x70c3, 0x70c4, 0x70c5, 0x70c6, 0x70c7, 0x7100, 0x7101,
                        0x7102, 0x7103, 0x7104, 0x7105, 0x7106, 0x7107, 0x7140, 0x7141, 0x7142, 0x7143, 0x7144, 0x7145,
                        0x7146, 0x7147, 0x7180, 0x7181, 0x7182, 0x7183, 0x7184, 0x7185, 0x7186, 0x7187, 0x71c0, 0x71c1,
                        0x71c2, 0x71c3, 0x71c4, 0x71c5, 0x71c6, 0x71c7
                };

        static constexpr uint_fast16_t morton_3_bitwise_y_512[512] =
                {0x0, 0x8, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38, 0x200, 0x208, 0x210, 0x218, 0x220, 0x228, 0x230, 0x238,
                 0x400, 0x408, 0x410, 0x418, 0x420, 0x428, 0x430, 0x438, 0x600, 0x608, 0x610, 0x618, 0x620, 0x628,
                 0x630, 0x638, 0x800, 0x808, 0x810, 0x818, 0x820, 0x828, 0x830, 0x838, 0xa00, 0xa08, 0xa10, 0xa18,
                 0xa20, 0xa28, 0xa30, 0xa38, 0xc00, 0xc08, 0xc10, 0xc18, 0xc20, 0xc28, 0xc30, 0xc38, 0xe00, 0xe08,
                 0xe10, 0xe18, 0xe20, 0xe28, 0xe30, 0xe38, 0x8000, 0x8008, 0x8010, 0x8018, 0x8020, 0x8028, 0x8030,
                 0x8038, 0x8200, 0x8208, 0x8210, 0x8218, 0x8220, 0x8228, 0x8230, 0x8238, 0x8400, 0x8408, 0x8410, 0x8418,
                 0x8420, 0x8428, 0x8430, 0x8438, 0x8600, 0x8608, 0x8610, 0x8618, 0x8620, 0x8628, 0x8630, 0x8638, 0x8800,
                 0x8808, 0x8810, 0x8818, 0x8820, 0x8828, 0x8830, 0x8838, 0x8a00, 0x8a08, 0x8a10, 0x8a18, 0x8a20, 0x8a28,
                 0x8a30, 0x8a38, 0x8c00, 0x8c08, 0x8c10, 0x8c18, 0x8c20, 0x8c28, 0x8c30, 0x8c38, 0x8e00, 0x8e08, 0x8e10,
                 0x8e18, 0x8e20, 0x8e28, 0x8e30, 0x8e38, 0x10000, 0x10008, 0x10010, 0x10018, 0x10020, 0x10028, 0x10030,
                 0x10038, 0x10200, 0x10208, 0x10210, 0x10218, 0x10220, 0x10228, 0x10230, 0x10238, 0x10400, 0x10408,
                 0x10410, 0x10418, 0x10420, 0x10428, 0x10430, 0x10438, 0x10600, 0x10608, 0x10610, 0x10618, 0x10620,
                 0x10628, 0x10630, 0x10638, 0x10800, 0x10808, 0x10810, 0x10818, 0x10820, 0x10828, 0x10830, 0x10838,
                 0x10a00, 0x10a08, 0x10a10, 0x10a18, 0x10a20, 0x10a28, 0x10a30, 0x10a38, 0x10c00, 0x10c08, 0x10c10,
                 0x10c18, 0x10c20, 0x10c28, 0x10c30, 0x10c38, 0x10e00, 0x10e08, 0x10e10, 0x10e18, 0x10e20, 0x10e28,
                 0x10e30, 0x10e38, 0x18000, 0x18008, 0x18010, 0x18018, 0x18020, 0x18028, 0x18030, 0x18038, 0x18200,
                 0x18208, 0x18210, 0x18218, 0x18220, 0x18228, 0x18230, 0x18238, 0x18400, 0x18408, 0x18410, 0x18418,
                 0x18420, 0x18428, 0x18430, 0x18438, 0x18600, 0x18608, 0x18610, 0x18618, 0x18620, 0x18628, 0x18630,
                 0x18638, 0x18800, 0x18808, 0x18810, 0x18818, 0x18820, 0x18828, 0x18830, 0x18838, 0x18a00, 0x18a08,
                 0x18a10, 0x18a18, 0x18a20, 0x18a28, 0x18a30, 0x18a38, 0x18c00, 0x18c08, 0x18c10, 0x18c18, 0x18c20,
                 0x18c28, 0x18c30, 0x18c38, 0x18e00, 0x18e08, 0x18e10, 0x18e18, 0x18e20, 0x18e28, 0x18e30, 0x18e38,
                 0x20000, 0x20008, 0x20010, 0x20018, 0x20020, 0x20028, 0x20030, 0x20038, 0x20200, 0x20208, 0x20210,
                 0x20218, 0x20220, 0x20228, 0x20230, 0x20238, 0x20400, 0x20408, 0x20410, 0x20418, 0x20420, 0x20428,
                 0x20430, 0x20438, 0x20600, 0x20608, 0x20610, 0x20618, 0x20620, 0x20628, 0x20630, 0x20638, 0x20800,
                 0x20808, 0x20810, 0x20818, 0x20820, 0x20828, 0x20830, 0x20838, 0x20a00, 0x20a08, 0x20a10, 0x20a18,
                 0x20a20, 0x20a28, 0x20a30, 0x20a38, 0x20c00, 0x20c08, 0x20c10, 0x20c18, 0x20c20, 0x20c28, 0x20c30,
                 0x20c38, 0x20e00, 0x20e08, 0x20e10, 0x20e18, 0x20e20, 0x20e28, 0x20e30, 0x20e38, 0x28000, 0x28008,
                 0x28010, 0x28018, 0x28020, 0x28028, 0x28030, 0x28038, 0x28200, 0x28208, 0x28210, 0x28218, 0x28220,
                 0x28228, 0x28230, 0x28238, 0x28400, 0x28408, 0x28410, 0x28418, 0x28420, 0x28428, 0x28430, 0x28438,
                 0x28600, 0x28608, 0x28610, 0x28618, 0x28620, 0x28628, 0x28630, 0x28638, 0x28800, 0x28808, 0x28810,
                 0x28818, 0x28820, 0x28828, 0x28830, 0x28838, 0x28a00, 0x28a08, 0x28a10, 0x28a18, 0x28a20, 0x28a28,
                 0x28a30, 0x28a38, 0x28c00, 0x28c08, 0x28c10, 0x28c18, 0x28c20, 0x28c28, 0x28c30, 0x28c38, 0x28e00,
                 0x28e08, 0x28e10, 0x28e18, 0x28e20, 0x28e28, 0x28e30, 0x28e38, 0x30000, 0x30008, 0x30010, 0x30018,
                 0x30020, 0x30028, 0x30030, 0x30038, 0x30200, 0x30208, 0x30210, 0x30218, 0x30220, 0x30228, 0x30230,
                 0x30238, 0x30400, 0x30408, 0x30410, 0x30418, 0x30420, 0x30428, 0x30430, 0x30438, 0x30600, 0x30608,
                 0x30610, 0x30618, 0x30620, 0x30628, 0x30630, 0x30638, 0x30800, 0x30808, 0x30810, 0x30818, 0x30820,
                 0x30828, 0x30830, 0x30838, 0x30a00, 0x30a08, 0x30a10, 0x30a18, 0x30a20, 0x30a28, 0x30a30, 0x30a38,
                 0x30c00, 0x30c08, 0x30c10, 0x30c18, 0x30c20, 0x30c28, 0x30c30, 0x30c38, 0x30e00, 0x30e08, 0x30e10,
                 0x30e18, 0x30e20, 0x30e28, 0x30e30, 0x30e38, 0x38000, 0x38008, 0x38010, 0x38018, 0x38020, 0x38028,
                 0x38030, 0x38038, 0x38200, 0x38208, 0x38210, 0x38218, 0x38220, 0x38228, 0x38230, 0x38238, 0x38400,
                 0x38408, 0x38410, 0x38418, 0x38420, 0x38428, 0x38430, 0x38438, 0x38600, 0x38608, 0x38610, 0x38618,
                 0x38620, 0x38628, 0x38630, 0x38638, 0x38800, 0x38808, 0x38810, 0x38818, 0x38820, 0x38828, 0x38830,
                 0x38838, 0x38a00, 0x38a08, 0x38a10, 0x38a18, 0x38a20, 0x38a28, 0x38a30, 0x38a38, 0x38c00, 0x38c08,
                 0x38c10, 0x38c18, 0x38c20, 0x38c28, 0x38c30, 0x38c38, 0x38e00, 0x38e08, 0x38e10, 0x38e18, 0x38e20,
                 0x38e28, 0x38e30, 0x38e38};


        template<typename morton, typename coord>
        static inline morton interleave_bitwise(const coord x, const coord y) {
            morton answer = 0;
            const static morton EIGHTBITMASK = 0x000000FF;
            for (unsigned int i = sizeof(coord); i > 0; --i) {
                unsigned int shift = (i - 1) * 8;
                answer =
                        answer << 16 |
                        morton_y_256[(y >> shift) & EIGHTBITMASK] |
                        morton_x_256[(x >> shift) & EIGHTBITMASK];
            }
            return answer;
        }

        template<typename morton, typename coord>
        static inline morton interleave_2_bitwise(const coord x, const coord y) {
            morton answer = 0;
            const static morton EIGHTBITMASK = 0x000000FF;
            for (unsigned int i = sizeof(coord); i > 0; --i) {
                unsigned int shift = (i - 1) * 8;
                answer =
                        answer << 16 |
                        morton_2_bitwise_y_256[(y >> shift) & EIGHTBITMASK] |
                        morton_2_bitwise_x_256[(x >> shift) & EIGHTBITMASK];
            }
            return answer;
        }

        //test for 64 Bit
        template<typename morton, typename coord>
        static inline morton interleave_3_bitwise(const coord x, const coord y) {
            morton answer = 0;
            const static morton NINEBITMASK = 0x000001FF;
            for (unsigned int i = sizeof(coord); i > 0; --i) {
                unsigned int shift = (i - 1) * 9;
                answer =
                        answer << 18 |
                        morton_3_bitwise_y_512[(y >> shift) & NINEBITMASK] |
                        morton_3_bitwise_x_512[(x >> shift) & NINEBITMASK];
            }
            return answer;
        }


        template<uint8_t t_k>
        struct interleave{
            static bool bits(const std::pair<uint, uint> lhs, const std::pair<uint, uint> rhs){
                throw std::runtime_error("not yet implemented");
            }
        };

        template<>
        struct interleave<2>{
            static bool inline bits(const std::pair<uint32_t, uint32_t> lhs, const std::pair<uint32_t, uint32_t> rhs){
                return interleave_bitwise<uint64_t, uint32_t>(lhs.second, lhs.first) < interleave_bitwise<uint64_t, uint32_t>(rhs.second, rhs.first);
            }
        };

        template<>
        struct interleave<4>{
            static bool inline bits(const std::pair<uint32_t, uint32_t> lhs, const std::pair<uint32_t, uint32_t> rhs){
                return interleave_2_bitwise<uint64_t, uint32_t>(lhs.second, lhs.first) < interleave_2_bitwise<uint64_t, uint32_t>(rhs.second, rhs.first);
            }
        };

        template<>
        struct interleave<8>{
            static bool inline bits(const std::pair<uint32_t, uint32_t> lhs, const std::pair<uint32_t, uint32_t> rhs){
                return interleave_3_bitwise<uint64_t, uint32_t>(lhs.second, lhs.first) < interleave_3_bitwise<uint64_t, uint32_t>(rhs.second, rhs.first);
            }
        };

        //corresponding_subtreeulate corresponding subtree on given level efficiently
        template<uint8_t t_k>
        struct access_shortcut_helper {
            static uint corresponding_subtree(uint32_t, uint32_t, uint8_t, uint) {
                //FIXME: add generic implementation/think about it
                throw new std::runtime_error("not yet implemented");
            }

            static uint64_t corresponding_subtree(uint64_t, uint64_t, uint8_t, uint) {
                //FIXME: add generic implementation/think about it
                throw new std::runtime_error("not yet implemented");
            }
        };

        template<>
        struct access_shortcut_helper<2> {
            static uint corresponding_subtree(uint32_t p, uint32_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 16) {
                    return interleaveFirstBits(p, q, level, m_real_size_of_max_element);
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit), look at https://github.com/Forceflow/libmorton
                    throw new std::runtime_error("not yet implemented");
                }
            }

            static uint64_t
            corresponding_subtree(uint64_t p, uint64_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 16) {
                    return interleaveFirstBits(p, q, level, m_real_size_of_max_element);
                } else if (level < 32) {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            /**
            * Interleaves the top h bits starting from bit l
            */
            static inline uint64_t interleaveFirstBits(uint64_t x, uint64_t y, int h, int l) {
                x = x >> (l - h);
                y = y >> (l - h);

                uint64_t z = ((morton_y_256[y >> 8] | morton_x_256[x >> 8]) << 16) |
                             morton_y_256[y & 0xFF] |
                             morton_x_256[x & 0xFF];
                return z;
            }
        };


        template<>
        struct access_shortcut_helper<4> {

            static uint corresponding_subtree(uint32_t p, uint32_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 8) {
                    return interleaveFirstBits(p, q, level, m_real_size_of_max_element);
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            static uint64_t
            corresponding_subtree(uint64_t p, uint64_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 8) {
                    return (interleaveFirstBits(p, q, level, m_real_size_of_max_element));
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            /**
            * Interleaves the top h bits starting from bit l pairwise
            */
            static inline uint32_t interleaveFirstBits(uint64_t x, uint64_t y, int h, int l) {
                x = x >> (l - (h * 2));
                y = y >> (l - (h * 2));

                uint32_t z = (morton_2_bitwise_y_256[y >> 8] | morton_2_bitwise_x_256[x >> 8]) << 16 |
                             morton_2_bitwise_y_256[y & 0xFF] |
                             morton_2_bitwise_x_256[x & 0xFF];
                return z;
            }
        };

        template<>
        struct access_shortcut_helper<8> {

            static uint corresponding_subtree(uint32_t p, uint32_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 5) {//maximum 18 Bit per coord for current level --> 36 Bit in Total
                    return interleaveFirstBits(p, q, level, m_real_size_of_max_element);
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            static uint64_t
            corresponding_subtree(uint64_t p, uint64_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 6) {
                    return (interleaveFirstBits(p, q, level, m_real_size_of_max_element));
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            /**
            * Interleaves the top h bits starting from bit l pairwise
            */
            static inline uint64_t interleaveFirstBits(uint64_t x, uint64_t y, int h, int l) {
                x = x >> (l - (h * 3));
                y = y >> (l - (h * 3));

                //only look at top 7 bits and later at bottom 9 bits, table gives 15 bits, which should be 18 bits (3 zeros at the end missing), shift the bits 3 to the right to obtain the 18 bits, then shift the 18 bits to the right in the case of 32 Bit discarding the top two bits (3+16 = 19)
                uint64_t z = (morton_3_bitwise_y_512[y >> 9] | morton_3_bitwise_x_512[x >> 9])  << 18;

                //std::cout << std::bitset<36>(z) << std::endl;


                //std::cout << std::bitset<36>(z) << std::endl;
                z = z | morton_3_bitwise_y_512[y & 0x1FF];

                //std::cout << std::bitset<36>(z) << std::endl;
                z = z | morton_3_bitwise_x_512[x & 0x1FF];

                //std::cout << std::bitset<36>(z) << std::endl;
                return z;
            }
        };

        template<>
        struct access_shortcut_helper<16> {

            static uint corresponding_subtree(uint32_t p, uint32_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 4) {
                    return interleaveFirstBits(p, q, level, m_real_size_of_max_element);
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            static uint64_t
            corresponding_subtree(uint64_t p, uint64_t q, uint8_t m_real_size_of_max_element, uint level) {
                if (level <= 4) {
                    return interleaveFirstBits(p, q, level, m_real_size_of_max_element);
                } else {
                    //FIXME: interleave top level+1 (more than 16 Bit)
                    throw new std::runtime_error("not yet implemented");
                }
            }

            /**
           * Interleaves the top h bits starting from bit l pairwise
           */
            static inline uint64_t interleaveFirstBits(uint64_t x, uint64_t y, int h, int l) {
                x = x >> (l - (h * 4)); //per Level 4 Bits per Coordinate (2^4 = 16)
                y = y >> (l - (h * 4));

                uint64_t z = ((y << 16) & 0xF0000000) |
                             ((x << 12) & 0x0F000000) |
                             ((y << 12) & 0x00F00000) |
                             ((x << 8) & 0x000F0000) |
                             ((y << 8) & 0x0000F000) |
                             ((x << 4) & 0x00000F00) |
                             ((y << 4) & 0x000000F0) |
                             (x & 0x0000000F);
                return z;
            }
        };

        //Efficient exponentiation by https://gist.github.com/orlp/3551590
        int64_t ipow(int64_t base, uint8_t exp) {
            static const uint8_t highest_bit_set[] = {
                    0, 1, 2, 2, 3, 3, 3, 3,
                    4, 4, 4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5, 5, 5, 5,
                    5, 5, 5, 5, 5, 5, 5, 5,
                    6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 6, 6, 6, 6, 255, // anything past 63 is a guaranteed overflow with base > 1
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
                    255, 255, 255, 255, 255, 255, 255, 255,
            };

            uint64_t result = 1;

            switch (highest_bit_set[exp]) {
                case 255: // we use 255 as an overflow marker and return 0 on overflow/underflow
                    if (base == 1) {
                        return 1;
                    }

                    if (base == -1) {
                        return 1 - 2 * (exp & 1);
                    }

                    return 0;
                case 6:
                    if (exp & 1) result *= base;
                    exp >>= 1;
                    base *= base;
                case 5:
                    if (exp & 1) result *= base;
                    exp >>= 1;
                    base *= base;
                case 4:
                    if (exp & 1) result *= base;
                    exp >>= 1;
                    base *= base;
                case 3:
                    if (exp & 1) result *= base;
                    exp >>= 1;
                    base *= base;
                case 2:
                    if (exp & 1) result *= base;
                    exp >>= 1;
                    base *= base;
                case 1:
                    if (exp & 1) result *= base;
                default:
                    return result;
            }
        }

        bool has_ending(std::string const &fullString, std::string const &ending) {
            if (fullString.length() >= ending.length()) {
                return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
            } else {
                return false;
            }
        }

        template<typename t_tv>
        uint64_t get_maximum(const t_tv &v) {
            using namespace k2_treap_ns;
            if (v.size() == 0) {
                return 0;
            }

            using t_e = typename t_tv::value_type;
            auto tupmax = [](t_e a) {
                return std::max(a.first, a.second);
            };
            auto max_it = std::max_element(std::begin(v), std::end(v), [&](t_e a, t_e b) {
                return tupmax(a) < tupmax(b);
            });
            return tupmax(*max_it);

        }

    } // end namepsace k2_treap_ns

} // end nomespace sdsl
#endif
