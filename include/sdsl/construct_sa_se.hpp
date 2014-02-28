#ifndef SDSL_CONSTRUCT_SA_SE
#define SDSL_CONSTRUCT_SA_SE

#include "io.hpp"
#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace sdsl
{

template<class int_vector_type>
uint64_t _get_next_lms_position(int_vector_type& text, uint64_t i)
{
    if (i >= text.size()-3) {
        return text.size()-1;
    }
    // text[i] is S-TYP or L-TYP
    uint64_t ci=text[i], cip1=text[i+1];
    while (ci <= cip1) {
        ++i;
        ci = cip1;
        cip1 = text[i+1];
    }
    // text[i] is L-TYP or S-TYP
    uint64_t candidate = i+1;
    while (ci >= cip1) {
        if (ci > cip1) {
            if (i+1==text.size()-1) {
                return text.size()-1;
            }
            candidate = i+1;
        }
        ++i;
        ci = cip1;
        cip1 = text[i+1];
    }
    return candidate;
}

void
_construct_sa_IS(int_vector<> &text, int_vector<> &sa,
                 std::string& filename_sa, size_t n, size_t text_offset,
                 size_t sigma, uint64_t recursion);

template <class int_vector_type >
void _construct_sa_se(int_vector_type& text, std::string filename_sa, uint64_t sigma, uint64_t recursion)
{
    std::string filename_text = tmp_file(filename_sa, "_text_rec"+util::to_string(recursion));
    store_to_file(text, filename_text);
    uint64_t n = text.size();
    uint64_t nsize = bits::hi(n)+1;
    uint8_t int_width = bits::hi(n-1)+1;
    uint64_t buffersize = 1024*1024/8;

    // Step 1 - Scan Text from right to left and count LMS, S and L characters and store lms_positions

    // Define variables
    size_t first_lms_pos=0;
    size_t number_of_lms_strings = 0;
    size_t bkt_s_last = 0, bkt_s_sum=0, bound_s=0, bkt_l_sum=0;
    int_vector<> C(sigma, 0, int_width);
    int_vector<> bkt_lms(sigma, 0, int_width);
    int_vector<> bkt_s(sigma, 0, int_width);
    int_vector<> bkt_l(sigma, 0, int_width);
    std::string filename_lms_pos_b = tmp_file(filename_sa, "_lms_pos_b"+util::to_string(recursion));
    size_t parts = 10;

    {
        int_vector_buffer<1> lms_pos_b(filename_lms_pos_b, std::ios::out, buffersize, 1);
        uint64_t ci = text[n-1];
        ++C[ci];
        bool was_s_typ = 1;
        for (size_t i=n-2; i<n; --i) {
            uint64_t cip1 = ci;
            ci = text[i];
            ++C[ci];
            if (was_s_typ) {
                ++bkt_s[text[i+1]];
                if (ci>cip1) {
                    ++bkt_lms[cip1];
                    lms_pos_b[i+1] = 1;
                    ++number_of_lms_strings;
                    first_lms_pos = i+1;
                    was_s_typ = 0;
                }
            } else if (ci<cip1) {
                was_s_typ = 1;
            }
        }
        if (was_s_typ) {
            ++bkt_s[ci];
        }
        bkt_l[0] = C[0]-bkt_s[0];
        for (size_t i=1; i<C.size(); ++i) {
            bkt_l[i] = C[i]-bkt_s[i];
            C[i] = C[i]+C[i-1];
        }
        lms_pos_b.close();
    }

    // Step 2 - Scan Text from right to left and detect LMS-Positions. Sort and write them to disk
    int_vector_buffer<> right(tmp_file(filename_sa, "_right"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
    size_t right_pointer=0;
    int_vector_buffer<> left(tmp_file(filename_sa, "_left"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
    size_t left_pointer=0;
    {
        for (size_t i=0, tmp2=0, tmp=0; i<sigma; ++i) {
            tmp += bkt_lms[i];
            bkt_lms[i] = tmp2;
            tmp2 = tmp;
        }
        int_vector_buffer<> lms_positions(tmp_file(filename_sa, "_lms_positions"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
        for (size_t i=n-2, was_s_typ=1, ci=text[n-1]; i<n; --i) {
            uint64_t cip1 = ci;
            ci = text[i];
            if (ci>cip1) {
                if (was_s_typ) {
                    lms_positions.push_back(bkt_lms[cip1]);
                    lms_positions.push_back(i+1);
                    ++bkt_lms[cip1];
                    was_s_typ = 0;
                }
            } else if (ci<cip1) {
                was_s_typ = 1;
            }
        }
        util::clear(text);
        {
            // Order lms_positions according to first character
            int_vector<> lms_strings(number_of_lms_strings, 0, int_width);
            for (size_t i=0; i<lms_positions.size();) {
                size_t idx = lms_positions[i++];
                size_t val = lms_positions[i++];
                lms_strings[idx] = val;
            }
            lms_positions.close(true);
            // Store it to file
            left_pointer = 0;
            for (size_t i=0; i<number_of_lms_strings; ++i) {
                left[left_pointer++] = lms_strings[number_of_lms_strings-i-1];
            }
        }
        load_from_file(text, filename_text);
    }
    left_pointer--;

    // Step 3 - Scan virtual array from left to right
    {
        // prepare bkt_l and backup it into bkt_lms
        for (size_t i=0, tmp=0; i<sigma; ++i) {
            tmp = bkt_l[i];
            bkt_l[i] = bkt_l_sum;
            bkt_l_sum += tmp;
            bkt_lms[i] = bkt_l[i];
        }

        size_t partsize = bkt_l_sum/parts+1;

        int_vector<> array(partsize, 0, int_width);
        std::vector< int_vector_buffer<> > cached_array(parts-1);
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i] = int_vector_buffer<>(tmp_file(filename_sa, "_rightbuffer"+util::to_string(i)+"_"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
        }

        for (size_t c=0, pos=0, offset=0; c<sigma; ++c) {
            // begin with array
            for (; pos<bkt_l[c]; ++pos) {
                // Load lazy values
                if (pos-offset >= partsize) {
                    offset += partsize;
                    for (size_t i=0, cur_part=pos/partsize-1; i<cached_array[cur_part].size();) {
                        size_t src = cached_array[cur_part][i++];
                        size_t val = cached_array[cur_part][i++];
                        array[src-offset] = val;
                    }
                    cached_array[pos/partsize-1].reset();
                }

                size_t idx = array[pos-offset];
                if (idx == 0) {
                    right[right_pointer++] = idx;
                } else {
                    size_t symbol = text[idx-1];
                    if (symbol >= c) {
                        size_t val = idx-1;
                        size_t src = bkt_l[symbol];
                        bkt_l[symbol] = bkt_l[symbol] + 1;
                        if ((src-offset)/partsize == 0) {
                            array[src-offset] = val;
                        } else {
                            size_t part = src/partsize-1;
                            cached_array[part].push_back(src);
                            cached_array[part].push_back(val);
                        }
                    } else {
                        right[right_pointer++] = idx;
                    }
                }
            }
            // continue with stack
            while (left_pointer < number_of_lms_strings and text[left[left_pointer]] == c) {
                size_t idx = left[left_pointer--];
                --idx;
                size_t symbol = text[idx];

                size_t val = idx;
                size_t src = bkt_l[symbol];
                bkt_l[symbol] = bkt_l[symbol] + 1;
                if ((src-offset)/partsize == 0) {
                    array[src-offset] = val;
                } else {
                    size_t part = src/partsize-1;
                    cached_array[part].push_back(src);
                    cached_array[part].push_back(val);
                }
            }
        }
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i].close(true);
        }

        // Restore bkt_l from bkt_lms
        for (size_t i=0; i<sigma; ++i) {
            bkt_l[i] = bkt_lms[i];
        }
    }
    right_pointer--;

    // Step 4 - Scan virtual array from right to left
    left_pointer = 0;
    left.reset();
    {
        // Prepare bkt_s and backup it into bkt_lms
        bkt_s_last=0, bkt_s_sum=0;
        for (size_t i=0; i<sigma; ++i) {
            bkt_s_sum += bkt_s[i];
            if (bkt_s[i]) {
                bkt_s[i] = bkt_s_sum;
                bkt_s_last = bkt_s_sum;
            } else {
                bkt_s[i] = bkt_s_sum;
            }
            bkt_lms[i] = bkt_s[i];
        }
        bound_s = bkt_s_sum;

        // Determine splitting parameters
        for (size_t i=0; i<sigma; ++i) {
            if (bkt_s[i] > bkt_s_sum/2) {
                bkt_s_sum = bkt_s[i];
                break;
            }
        }

        size_t partsize = bound_s/parts+1;
        int_vector<> array(partsize, 0, int_width);
        std::vector< int_vector_buffer<> > cached_array(parts-1);
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i] = int_vector_buffer<>(tmp_file(filename_sa, "_leftbuffer"+util::to_string(i)+"_"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
        }
        for (size_t c=sigma-1, pos=bkt_s_last-1, offset=partsize*(parts-1); c<sigma; --c) {
            // begin with array
            for (; pos+1 > bkt_s[c]; --pos) {
                while (pos < offset) {
                    // Load buffered values
                    offset -= partsize;
                    for (size_t i=0, cur_part=offset/partsize; i<cached_array[cur_part].size();) {
                        size_t src = cached_array[cur_part][i++];
                        size_t val = cached_array[cur_part][i++];
                        array[src-offset] = val;
                    }
                    cached_array[offset/partsize].reset();
                }

                size_t idx = array[pos-offset];
                if (idx==0) {
                    idx = n;
                }
                --idx;
                size_t symbol = text[idx];
                if (symbol <= c) {
                    bkt_s[symbol] = bkt_s[symbol] - 1;
                    size_t val = idx;
                    size_t src = bkt_s[symbol];
                    if (src >= offset) {
                        array[src-offset] = val;
                    } else {
                        size_t part = src/partsize;
                        cached_array[part].push_back(src);
                        cached_array[part].push_back(val);
                    }
                } else {
                    left[left_pointer++] = array[pos-offset];
                }
            }

            // continue with stack
            while (right_pointer < number_of_lms_strings and text[right[right_pointer]] == c) {
                size_t idx = right[right_pointer--];
                if (idx == 0) {
                    idx = n;
                }
                --idx;
                size_t symbol = text[idx];
                bkt_s[symbol] = bkt_s[symbol] - 1;

                size_t val = idx;
                size_t src = bkt_s[symbol];
                if (src >= offset) {
                    array[src-offset] = val;
                } else {
                    size_t part = src/partsize;
                    cached_array[part].push_back(src);
                    cached_array[part].push_back(val);
                }
            }
        }
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i].close(true);
        }
        // Restore bkt_s from bkt_lms
        for (size_t i=0; i<sigma; ++i) {
            bkt_s[i] = bkt_lms[i];
        }
    }
    right.buffersize(0);
    right.reset();
    right_pointer = 0;
    --left_pointer;

    // Step 5 - Detect same lms-Strings, write text to file
    int_vector<1> same_lms(number_of_lms_strings, false);
    size_t last_end_pos = first_lms_pos, order = number_of_lms_strings-1;
    same_lms[number_of_lms_strings-1] = true;
    for (size_t i=number_of_lms_strings-2, a=0, b=0, last_a=left[number_of_lms_strings-1]; i<number_of_lms_strings; --i) {
        b = last_a;
        a = left[i];
        last_a = a;

        size_t end_pos = _get_next_lms_position(text, a);
        if (end_pos-a == last_end_pos-b) {
            while (a < end_pos and text[a] == text[b]) {
                ++a;
                ++b;
            }
            if (text[a] == text[b]) {
                same_lms[i] = true;
                --order;
            }
        }
        last_end_pos = end_pos;
    }
    util::clear(text);

    // Step 7 - Create renamed string
    int_vector<> text_rec;
    if (recursion==0) {
        text_rec.width((bits::hi(order+1)+1));
    } else {
        text_rec.width((bits::hi(number_of_lms_strings+1)+1));
    }
    text_rec.resize(number_of_lms_strings);
    util::_set_zero_bits(text_rec);
    {
        if (recursion==0 and n/2*text_rec.width()>8*n) {
            size_t size_of_part = n/4+3;
            text_rec.resize(size_of_part);
            util::_set_zero_bits(text_rec);
            order = 0;
            for (size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i) {
                if (!same_lms[i]) {
                    ++order;
                }
                if (left[i]/2 >= size_of_part) {
                    text_rec[(left[i]/2)-size_of_part] = order;
                }
            }
            std::string filename_text_rec_part2 = tmp_file(filename_sa, "_text_rec_part2"+util::to_string(recursion));
            size_t pos = 0;
            for (size_t i=0; i<size_of_part; ++i) {
                if (text_rec[i]>0) {
                    text_rec[pos++] = text_rec[i];
                }
            }
            text_rec.resize(pos);
            store_to_file(text_rec, filename_text_rec_part2);
            text_rec.resize(size_of_part);
            util::_set_zero_bits(text_rec);
            order = 0;
            for (size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i) {
                if (!same_lms[i]) {
                    ++order;
                }
                if (left[i]/2 < size_of_part) {
                    text_rec[left[i]/2] = order;
                }
            }
            pos = 0;
            for (size_t i=0; i<size_of_part; ++i) {
                if (text_rec[i]>0) {
                    text_rec[pos++] = text_rec[i];
                }
            }
            text_rec.resize(number_of_lms_strings);
            int_vector_buffer<> buf(filename_text_rec_part2, std::ios::in, 1024*1024);
            for (size_t i=0; i<buf.size(); ++i) {
                text_rec[pos++] = buf[i];
            }
            buf.close(true);
            text_rec[number_of_lms_strings-1] = 0;
        } else {
            text_rec.resize(n/2+1);
            util::_set_zero_bits(text_rec);
            order = 0;
            for (size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i) {
                if (!same_lms[i]) {
                    ++order;
                }
                text_rec[left[left_pointer--]/2] = order;
            }
            for (size_t i=0, pos=0; i<text_rec.size(); ++i) {
                if (text_rec[i]>0) {
                    text_rec[pos++] = text_rec[i];
                }
            }
            text_rec[number_of_lms_strings-1] = 0;
            text_rec.resize(number_of_lms_strings);
        }
    }
    util::clear(same_lms);
    left.buffersize(0);
    left.reset();

    // Step 8 - Determine complete LMS order (recursivly)
    int_vector<> isa_rec;
    std::string filename_sa_rec = tmp_file(filename_sa, "_sa_rec"+util::to_string(recursion+1));
    if (text_rec.size() > order+1) {
        if (recursion==0) {
            memory_monitor::event("begin _construct_sa");
            _construct_sa_se<int_vector<> >(text_rec, filename_sa_rec, order+1, recursion+1);
            memory_monitor::event("end   _construct_sa");
        } else {
            text_rec.resize(text_rec.size()*2);
            for (size_t i=0; i<number_of_lms_strings; ++i) {
                text_rec[number_of_lms_strings+i] = text_rec[i];
                text_rec[i] = 0;
            }
            memory_monitor::event("begin sa_simple");
            _construct_sa_IS(text_rec, text_rec, filename_sa_rec, number_of_lms_strings, number_of_lms_strings, order+1, recursion+1);
            memory_monitor::event("end   sa_simple");
            // SA' in first half, S' in second half
            text_rec.resize(number_of_lms_strings);
            store_to_file(text_rec, filename_sa_rec);
        }
    } else {
        isa_rec.swap(text_rec);
    }
    util::clear(text_rec);

    // Step 9 - Initiate left for scan in step 12
    if (isa_rec.size() > 0) {
        // isa_rec exists //TODO test if its better to create sa_rec
        // TODO always enough memory? in memory: isa_rec, lms_pos_b, select_support, tmp_left, leftbuffer
        // load bit_vector lms_positions and create select support
        bit_vector lms_pos_b(n);
        load_from_file(lms_pos_b, filename_lms_pos_b);
        sdsl::remove(filename_lms_pos_b);
        select_support_mcl<> lms_select_support;                 // select_support for bit_vector
        util::init_support(lms_select_support, &lms_pos_b);  // Create select_support
        // write left
        int_vector<> tmp_left(number_of_lms_strings, 0, int_width);
        for (size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i) {
            size_t idx = isa_rec[i];
            size_t val = lms_select_support.select(i+1); //TODO test alternative without select support: look for 1 in lms_pos_b (backwards)
            tmp_left[idx] = val;
        }
        util::clear(lms_select_support);
        util::clear(lms_pos_b);
        util::clear(isa_rec);
        // write to left
        left.buffersize(buffersize);
        left_pointer = 0;
        for (; left_pointer<number_of_lms_strings; ++left_pointer) {
            left[left_pointer] = tmp_left[number_of_lms_strings-left_pointer-1];
        }
        left_pointer--;
        util::clear(tmp_left);
    } else {
        left.buffersize(buffersize);
        left_pointer = 0;
        {
            // load bit_vector lms_positions and create select support
            bit_vector lms_pos_b(n);
            load_from_file(lms_pos_b, filename_lms_pos_b);
            sdsl::remove(filename_lms_pos_b);
            select_support_mcl<> lms_select_support;                 // select_support for bit_vector
            util::init_support(lms_select_support, &lms_pos_b);      // create select_support
            // write to left sa_rec buffered
            int_vector_buffer<> sa_rec_buf(filename_sa_rec, std::ios::in, buffersize, nsize);
            for (uint64_t i=0; i<sa_rec_buf.size(); ++i) {
                uint64_t pos = lms_select_support.select(sa_rec_buf[i]+1);
                left[number_of_lms_strings-1-left_pointer++] = pos;
            }
            sa_rec_buf.close(true);
            left_pointer--;
        }
        //TODO test sa_rec unbuffered in recursion level 1 -> space still good?
    }


    // Step 10 - Reload text
    load_from_file(text, filename_text);
    sdsl::remove(filename_text);

    // Step 12 - Scan virtual array from left to right second time
    right.buffersize(buffersize);
    right_pointer = 0;
    int_vector_buffer<> cached_sa(filename_sa, std::ios::out, buffersize, nsize);
    size_t sa_pointer = 0;
    {
        size_t partsize = bkt_l_sum/parts+1;
        int_vector<> array(partsize, 0, int_width);
        std::vector< int_vector_buffer<> > cached_array(parts-1);
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i] = int_vector_buffer<>(tmp_file(filename_sa, "_rightbuffer"+util::to_string(i)+"_"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
        }

        for (size_t c=0, pos=0, offset=0; c<sigma; ++c) {
            // begin with array
            for (; pos<bkt_l[c]; ++pos) {
                // Load lazy values
                if (pos-offset >= partsize) {
                    offset += partsize;
                    for (size_t i=0, cur_part=pos/partsize-1; i<cached_array[cur_part].size();) {
                        size_t src = cached_array[cur_part][i++];
                        size_t val = cached_array[cur_part][i++];
                        array[src-offset] = val;
                    }
                    cached_array[pos/partsize-1].reset();
                }
                size_t idx = array[pos-offset];
                if (idx == 0) {
                    cached_sa[sa_pointer++] = idx;
                    right[right_pointer++] = idx;
                } else {
                    size_t symbol = text[idx-1];
                    cached_sa[sa_pointer++] = idx;
                    if (symbol >= c) {
                        size_t val = idx-1;
                        size_t src = bkt_l[symbol];
                        bkt_l[symbol] = bkt_l[symbol] + 1;
                        if ((src-offset)/partsize == 0) {
                            array[src-offset] = val;
                        } else {
                            size_t part = src/partsize-1;
                            cached_array[part].push_back(src);
                            cached_array[part].push_back(val);
                        }
                    } else {
                        right[right_pointer++] = idx;
                    }
                }
            }
            sa_pointer = C[c];
            // continue with stack
            while (left_pointer < number_of_lms_strings and text[left[left_pointer]] == c) {
                size_t idx = left[left_pointer--];
                if (idx == 0) {
                    idx = n;
                }
                --idx;
                size_t symbol = text[idx];
                size_t val = idx;
                size_t src = bkt_l[symbol];
                bkt_l[symbol] = bkt_l[symbol] + 1;
                if ((src-offset)/partsize == 0) {
                    array[src-offset] = val;
                } else {
                    size_t part = src/partsize-1;
                    cached_array[part].push_back(src);
                    cached_array[part].push_back(val);
                }
            }
        }
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i].close(true);
        }
    }
    left.close(true);
    right_pointer--;

    // Step 13 - Scan virtual array from right to left second time
    {
        size_t partsize = bound_s/parts+1;

        int_vector<> array(partsize, 0, int_width);
        std::vector< int_vector_buffer<> > cached_array(parts-1);
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i] = int_vector_buffer<>(tmp_file(filename_sa, "_leftbuffer"+util::to_string(i)+"_"+util::to_string(recursion)), std::ios::out, buffersize, nsize);
            //		cached_array_pointer[i] = 0;
        }
        for (size_t c=sigma-1, pos=bkt_s_last-1, offset=partsize*(parts-1); c<sigma; --c) {
            // begin with array
            assert(c < C.size());
            sa_pointer = C[c]-1;
            for (; pos+1 > bkt_s[c]; --pos) {
                while (pos < offset) {
                    // Load buffered values
                    offset -= partsize;
                    for (size_t i=0, cur_part=offset/partsize; i<cached_array[cur_part].size();) {
                        size_t src = cached_array[cur_part][i++];
                        size_t val = cached_array[cur_part][i++];
                        assert((src-offset) < array.size());
                        array[src-offset] = val;
                    }
                    assert((offset/partsize) < parts-1);
                    cached_array[offset/partsize].reset();
                }

                assert((pos-offset) < array.size());
                size_t idx = array[pos-offset];
                if (idx==0) {
                    idx = n;
                }
                --idx;
                assert((idx) < text.size());
                size_t symbol = text[idx];
                if (symbol <= c) {
                    if (idx==n-1) {
                        cached_sa[sa_pointer--] = 0;
                    } else {
                        cached_sa[sa_pointer--] = idx+1;
                    }
                    assert((symbol) < bkt_s.size());
                    bkt_s[symbol] = bkt_s[symbol] - 1;
                    size_t val = idx;
                    size_t src = bkt_s[symbol];
                    if (src >= offset) {
                        assert((src-offset) < array.size());
                        array[src-offset] = val;
                    } else {
                        size_t part = src/partsize;
                        assert(part < parts-1);
                        cached_array[part].push_back(src);
                        cached_array[part].push_back(val);
                    }
                } else {
                    if (idx==n-1) {
                        cached_sa[sa_pointer--] = 0;
                    } else {
                        cached_sa[sa_pointer--] = idx+1;
                    }
                }
            }
            // continue with stack
            while (right_pointer < number_of_lms_strings and text[right[right_pointer]] == c) {
                size_t idx = right[right_pointer--];
                if (idx == 0) {
                    idx = n;
                }
                --idx;
                size_t symbol = text[idx];
                assert((symbol) < bkt_s.size());
                bkt_s[symbol] = bkt_s[symbol] - 1;

                size_t val = idx;
                size_t src = bkt_s[symbol];
                if (src >= offset) {
                    assert((src-offset) < array.size());
                    array[src-offset] = val;
                } else {
                    size_t part = src/partsize;
                    assert((part) < parts-1);
                    cached_array[part].push_back(src);
                    cached_array[part].push_back(val);
                }
            }
        }
        for (size_t i=0; i<cached_array.size(); ++i) {
            cached_array[i].close(true);
        }
    }
    right.close(true);
    cached_sa.close();

    return;
}

} // end namespace

#endif
