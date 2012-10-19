/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
#include "sdsl/csa_uncompressed.hpp"

namespace sdsl
{

csa_uncompressed::csa_uncompressed(tMSS& file_map, const std::string& dir, const std::string& id):char2comp(m_char2comp),
	comp2char(m_comp2char), C(m_C), sigma(m_sigma), psi(m_psi), bwt(m_bwt), sa_sample(m_sa), isa_sample(m_isa) {

	int_vector_file_buffer<8> text_buf(file_map["text"].c_str());
	int_vector_file_buffer<>  sa_buf(file_map["sa"].c_str());
	size_type n = text_buf.int_vector_size;
	algorithm::set_text<csa_uncompressed>(text_buf, n, m_C, m_char2comp, m_comp2char, m_sigma);
	util::assign(m_sa, sa_sample_type(sa_buf));
	algorithm::set_isa_samples<csa_uncompressed>(sa_buf, m_isa);
	m_psi = psi_type(this);
	m_bwt = bwt_type(this);
	write_R_output("csa", "store ISA","begin",1,0);
	if (!util::store_to_file(m_isa, (dir+"isa_"+id).c_str(), true)) {
		throw std::ios_base::failure("#csa_uncompressed: Cannot store ISA to file system!");
	} else {
		file_map["isa"] = dir+"isa_"+id;
	}
	write_R_output("csa", "store ISA","end",1,0);
}

csa_uncompressed::const_iterator csa_uncompressed::begin()const
{
    return const_iterator(this, 0);
}


csa_uncompressed::const_iterator csa_uncompressed::end()const
{
    return const_iterator(this, size());
}

void csa_uncompressed::copy(const csa_uncompressed& csa)
{
    m_sa = csa.m_sa;
    m_isa = csa.m_isa;
    m_sigma		 = csa.m_sigma;
    m_char2comp  = csa.m_char2comp;
    m_comp2char  = csa.m_comp2char;
    m_C = csa.m_C;
    m_psi = psi_type(this);
    m_bwt = bwt_type(this);
}

csa_uncompressed& csa_uncompressed::operator=(const csa_uncompressed& csa)
{
    if (this != &csa) {
        copy(csa);
    }
    return *this;
}


csa_uncompressed::size_type csa_uncompressed::serialize(std::ostream& out)const
{
    size_type written_bytes = 0;
    written_bytes += m_sa.serialize(out);
    written_bytes += m_isa.serialize(out);
    written_bytes += m_char2comp.serialize(out);
    written_bytes += m_comp2char.serialize(out);
    written_bytes += m_C.serialize(out);
    written_bytes += util::write_member(m_sigma, out);
    return written_bytes;
}

void csa_uncompressed::load(std::istream& in)
{
    m_sa.load(in);
    m_isa.load(in);
    m_char2comp.load(in);
    m_comp2char.load(in);
    m_C.load(in);
    util::read_member(m_sigma, in);
    m_psi = psi_type(this);
    m_bwt = bwt_type(this);
}

void csa_uncompressed::swap(csa_uncompressed& csa)
{
    if (this != &csa) {
        m_sa.swap(csa.m_sa);
        m_isa.swap(csa.m_isa);
        m_psi.swap(csa.m_psi);
        m_char2comp.swap(csa.m_char2comp);
        m_comp2char.swap(csa.m_comp2char);
        m_C.swap(csa.m_C);
        std::swap(m_sigma, csa.m_sigma);
        m_psi = psi_type(this);
        m_bwt = bwt_type(this);
        csa.m_psi = psi_type(&csa);
        csa.m_bwt = bwt_type(&csa);
    }
}

bool csa_uncompressed::operator==(const csa_uncompressed& csa)const
{
    for (uint16_t i=0; i<256; ++i)
        if (m_char2comp[i] != csa.m_char2comp[i] or m_comp2char[i] != csa.m_comp2char[i] or m_C[i] != csa.m_C[i])
            return false;
    return m_sa == csa.m_sa and m_isa == csa.m_isa and m_psi == csa.m_psi and m_C[256] == csa.m_C[256] and m_sigma == csa.m_sigma;
}

bool csa_uncompressed::operator!=(const csa_uncompressed& csa)const
{
    return !(*this == csa);
}

csa_uncompressed::size_type csa_uncompressed::rank_bwt(size_type i, const unsigned char c)const
{
    // TODO: special case if c == BWT[i-1] we can use LF to get a constant time answer
    unsigned char cc = m_char2comp[c];
    if (cc==0 and c!=0)  // character is not in the text => return 0
        return 0;
    // binary search the interval [C[cc]..C[cc+1]-1] for the result
    size_type lower_b = m_C[cc], upper_b = m_C[cc+1]; // lower_b inclusive, upper_b exclusive
    while (lower_b+1 < upper_b) {
        size_type mid = (lower_b+upper_b)/2;
        if (m_psi[mid] >= i)
            upper_b = mid;
        else
            lower_b = mid;
    }
    if (lower_b > m_C[cc])
        return lower_b - m_C[cc] + 1;
    else { // lower_b == m_C[cc]
        return m_psi[lower_b] < i;// 1 if m_psi[lower_b]<i, 0 otherwise
    }
}

csa_uncompressed::size_type csa_uncompressed::select_bwt(size_type i, const unsigned char c)const
{
    assert(i > 0);
    unsigned char cc = m_char2comp[c];
    if (cc==0 and c!=0)  // character is not in the text => return size()
        return size();
    assert(cc != 255);
    if (m_C[cc]+i-1 <  m_C[cc+1]) {
        return m_psi[m_C[cc]+i-1];
    } else
        return size();
}

} // end namespace sdsl
