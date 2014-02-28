#include "sdsl/construct_sa_se.hpp"

namespace sdsl
{

void _construct_sa_IS(int_vector<> &text, int_vector<> &sa, std::string& filename_sa, size_t n, size_t text_offset, size_t sigma, uint64_t recursion)
{
    uint64_t buffersize = 1024*1024/8;

    size_t name = 0;
    size_t number_of_lms_strings = 0;
    std::string filename_c_array = tmp_file(filename_sa, "_c_array"+util::to_string(recursion));
    // Phase 1
    {
        std::vector<uint64_t> bkt(sigma, 0);
        // Step 1 - Count characters into c array
        // TODO: better create this in higher recursion-level
        for (size_t i=0; i<n; ++i) {
            ++bkt[text[text_offset+i]];
        }

        // Step 1.5 save them into cached_external_array
        int_vector_buffer<> c_array(filename_c_array, std::ios::out, buffersize, 64);
        for (size_t c=0; c<sigma; ++c) {
            c_array[c] = bkt[c];
        }

        // Step 2 Calculate End-Pointer of Buckets
        bkt[0] = 0;
        for (size_t c=1; c<sigma; ++c) {
            bkt[c] = bkt[c-1]+bkt[c];
        }

        // Step 3 - Insert S*-positions into correct bucket of SA but not in correct order inside the buckets
        for (size_t i=n-2, was_s_typ = 1; i<n; --i) {
            if (text[text_offset+i]>text[text_offset+i+1]) {
                if (was_s_typ) {
                    sa[bkt[text[text_offset+i+1]]--] = i+1;
                    ++number_of_lms_strings;
                    was_s_typ = 0;
                }
            } else if (text[text_offset+i]<text[text_offset+i+1]) {
                was_s_typ = 1;
            }
        }

        // Step 4 - Calculate Begin-Pointer of Buckets
        bkt[0] = 0;
        for (size_t c=1; c<sigma; ++c) {
            bkt[c] = bkt[c-1]+c_array[c-1];
        }

        // Step 5 - Scan from Left-To-Right to induce L-Types
        for (size_t i=0; i<n; ++i) {
            if (sa[i] > 0 and text[text_offset+ sa[i] ] <= text[text_offset+ sa[i]-1 ]) { // faster than if(sa[i]>0 and bkt_beg[text[ sa[i]-1 ]] > i)
                sa[bkt[text[text_offset+ sa[i]-1 ]]++] = sa[i]-1;
                sa[i] = 0;
            }
        }

        // Step 6 - Scan from Right-To-Left to induce S-Types
        bkt[0] = 0;
        for (size_t c=1; c<sigma; ++c) {
            bkt[c] = bkt[c-1]+c_array[c];
        }
        c_array.close();
        c_array.buffersize(0);

        for (size_t i=n-1, endpointer=n; i<n; --i) {
            if (sa[i]>0) {
                if (text[text_offset+ sa[i]-1 ] <= text[text_offset+ sa[i] ]) { // faster than if(bkt_end[text[ sa[i]-1 ]] < i)
                    sa[bkt[text[text_offset+ sa[i]-1 ]]--] = sa[i]-1;
                } else {
                    sa[--endpointer] = sa[i];
                }
                sa[i] = 0;
            }
        }

        // Step 7 - Determine length of LMS-Strings
        for (size_t i=n-2, end=n-2, was_s_typ = 1; i<n; --i) {
            if (text[text_offset+i]>text[text_offset+i+1]) {
                if (was_s_typ) {
                    sa[(i+1)>>1] = end-i;
                    end = i+1;
                    was_s_typ = 0;
                }
            } else if (text[text_offset+i]<text[text_offset+i+1]) {
                was_s_typ = 1;
            }
        }

        // Step 8 - Rename
        for (size_t i=n-number_of_lms_strings+1, cur_pos=0, cur_len=0, last_pos=n-1, last_len=1; i<n; ++i) {
            cur_pos = sa[i];
            cur_len = sa[(cur_pos>>1)];
            if (cur_len == last_len) {
                size_t l = 0;
                while (l < cur_len and text[text_offset+cur_pos+l] == text[text_offset+last_pos+l]) {
                    ++l;
                }
                if (l >= cur_len) {
                    --name;
                }
            }
            sa[(cur_pos>>1)] = ++name;
            last_pos = cur_pos;
            last_len = cur_len;
        }
    }

    // Step 9 - Calculate SA of new string - in most cases recursive
    if (name+1 < number_of_lms_strings) {
        // Move Names to the end
        for (size_t i=0, t=n-number_of_lms_strings; i<(n>>1); ++i) {
            if (sa[i] > 0) {
                sa[t++] = sa[i];
                sa[i] = 0;
            }
        }
        sa[n-1] = 0;

        // Recursive call
        std::string filename_sa_rec = tmp_file(filename_sa, "_sa_rec"+util::to_string(recursion+1));
        _construct_sa_IS(sa, sa, filename_sa_rec, number_of_lms_strings, n-number_of_lms_strings, name+1, recursion+1);

        for (size_t i=n-2, endpointer = n-1, was_s_typ = 1; i<n; --i) {
            if (text[text_offset+i]>text[text_offset+i+1]) {
                if (was_s_typ) {
                    sa[endpointer--] = i+1;
                    was_s_typ = 0;
                }
            } else if (text[text_offset+i]<text[text_offset+i+1]) {
                was_s_typ = 1;
            }
        }

        // Sort S*-positions in correct order into SA
        for (size_t i=0; i<number_of_lms_strings; ++i) {
            size_t pos = sa[i];
            sa[i] = sa[n-number_of_lms_strings+pos];
            sa[n-number_of_lms_strings+pos] = 0;
        }
    } else {
        // Move s*-Positions to front
        sa[0] = n-1;
        for (size_t i=1; i<number_of_lms_strings; ++i) {
            sa[i] = sa[n-number_of_lms_strings+i];
            sa[n-number_of_lms_strings+i] = 0;
        }
        // Clear lex. names
        for (size_t i=number_of_lms_strings; i<(n>>1); ++i) {
            sa[i] = 0;
        }
    }

    // Phase 3
    {
        // Step 10 - Count characters into c array

        // Step 11 - Calculate End-Pointer of Buckets
        int_vector_buffer<> c_array(filename_c_array, std::ios::in, buffersize, 64);
        std::vector<uint64_t> bkt(sigma, 0);
        for (size_t c=1; c<sigma; ++c) {
            bkt[c] = bkt[c-1]+c_array[c];
        }

        // Step 12 - Move S*-positions in correct order into SA
        for (size_t i=number_of_lms_strings-1; i<n; --i) {
            size_t pos = sa[i];
            sa[i] = 0;
            sa[ bkt[text[text_offset+pos]]-- ] = pos;
        }

        // Step 13 - Calculate Begin-Pointer of Buckets
        bkt[0] = 0;
        for (size_t c=1; c<sigma; ++c) {
            bkt[c] = bkt[c-1]+c_array[c-1];
        }

        // Step 14 - Scan from Left-To-Right to induce L-Types
        for (size_t i=0; i<n; ++i) {
            if (sa[i] > 0 and text[text_offset+ sa[i] ] <= text[text_offset+ sa[i]-1 ]) { // faster than if(sa[i]>0 and bkt_beg[text[ sa[i]-1 ]] > i)
                sa[bkt[text[text_offset+ sa[i]-1 ]]++] = sa[i]-1;
            }
        }

        // Step 15 - Scan from Right-To-Left to induce S-Types
        bkt[0] = 0;
        for (size_t c=1; c<sigma; ++c) {
            bkt[c] = bkt[c-1]+c_array[c];
        }
        for (size_t i=n-1; i<n; --i) {
            if (sa[i] > 0 and text[text_offset+sa[i]-1] <= text[text_offset+sa[i]]) {
                sa[bkt[text[text_offset+ sa[i]-1 ]]--] = sa[i]-1;
            }
        }
        c_array.close(true);
    }
}

}
