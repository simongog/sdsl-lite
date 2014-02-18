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
#include "sdsl/csa_alphabet_strategy.hpp"

namespace sdsl
{

const char* key_trait<8>::KEY_BWT = conf::KEY_BWT;
const char* key_trait<8>::KEY_TEXT = conf::KEY_TEXT;

byte_alphabet::byte_alphabet(int_vector_buffer<8>& text_buf, int_vector_size_type len):
    char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
    m_sigma = 0;
    if (0 == len or 0 == text_buf.size())
        return;
    assert(len <= text_buf.size());
    // initialize vectors
    util::assign(m_C	    , int_vector<64>(257, 0));
    util::assign(m_char2comp, int_vector<8>(256,0));
    util::assign(m_comp2char, int_vector<8>(256,0));
    // count occurrences of each symbol
    for (size_type i=0; i < len; ++i) {
        ++m_C[text_buf[i]];
    }
    assert(1 == m_C[0]); // null-byte should occur exactly once
    m_sigma = 0;
    for (int i=0; i<256; ++i)
        if (m_C[i]) {
            m_char2comp[i] 	 	= m_sigma;
            m_comp2char[sigma]  = i;
            m_C[m_sigma]		= m_C[i];
            ++m_sigma;
        }
    m_comp2char.resize(m_sigma);
    m_C.resize(m_sigma+1);
    for (int i=(int)m_sigma; i > 0; --i) m_C[i] = m_C[i-1];
    m_C[0] = 0;
    for (int i=1; i <= (int)m_sigma; ++i) m_C[i] += m_C[i-1];
    assert(C[sigma]==len);
}


byte_alphabet::byte_alphabet(): char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
    m_sigma = 0;
}

void byte_alphabet::copy(const byte_alphabet& bas)
{
    m_char2comp = bas.m_char2comp;
    m_comp2char = bas.m_comp2char;
    m_C			= bas.m_C;
    m_sigma		= bas.m_sigma;
}

byte_alphabet::byte_alphabet(const byte_alphabet& bas): char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
    copy(bas);
}

byte_alphabet& byte_alphabet::operator=(const byte_alphabet& bas)
{
    if (this != &bas) {
        copy(bas);
    }
    return *this;
}

byte_alphabet& byte_alphabet::operator=(byte_alphabet&& bas)
{
    if (this != &bas) {
        m_char2comp = std::move(bas.m_char2comp);
        m_comp2char = std::move(bas.m_comp2char);
        m_C         = std::move(bas.m_C);
        m_sigma     = std::move(bas.m_sigma);
    }
    return *this;
}

void byte_alphabet::swap(byte_alphabet& bas)
{
    m_char2comp.swap(bas.m_char2comp);
    m_comp2char.swap(bas.m_comp2char);
    m_C.swap(bas.m_C);
    std::swap(m_sigma, bas.m_sigma);
}

byte_alphabet::size_type byte_alphabet::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_char2comp.serialize(out, child, "m_char2comp");
    written_bytes += m_comp2char.serialize(out, child, "m_comp2char");
    written_bytes += m_C.serialize(out, child, "m_C");
    written_bytes += write_member(m_sigma, out, child, "m_sigma");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void byte_alphabet::load(std::istream& in)
{
    m_char2comp.load(in);
    m_comp2char.load(in);
    m_C.load(in);
    read_member(m_sigma, in);
}
}
