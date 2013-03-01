#include "sdsl/wt_helper.hpp"

namespace sdsl{

void calculate_character_occurences(int_vector_file_buffer<8>& text, const int_vector_size_type size, int_vector_size_type* C)
{
    text.reset();
    if (text.int_vector_size < size) {
        throw std::logic_error("calculate_character_occurrences: stream size is smaller than size!");
        return;
    }
    for (int_vector_size_type i=0, r_sum=0, r = text.load_next_block(); r_sum < size;) {
        if (r_sum + r > size) {  // read not more than size chars in the next loop
            r = size-r_sum;
        }
        for (; i < r_sum+r; ++i) {
            ++C[text[i-r_sum]];
        }
        r_sum += r; r = text.load_next_block();
    }
}


}
