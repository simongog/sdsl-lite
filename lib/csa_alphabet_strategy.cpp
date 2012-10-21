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

namespace sdsl{

	byte_alphabet_stategy::byte_alphabet_stategy(): char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
	{};


	void byte_alphabet_stategy::swap(byte_alphabet_stategy& bas){
		m_char2comp.swap(bas.m_char2comp);
		m_comp2char.swap(bas.m_comp2char);
		m_C.swap(bas.m_C);
		std::swap(m_sigma, bas.m_sigma);
	}

	typename byte_alphabet_stategy::size_type byte_alphabet_stategy::serialize(std::ostream& out, structure_tree_node*v, std::string name)const{
		structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		size_type written_bytes = 0;
        written_bytes += m_char2comp.serialize(out, child, "m_char2comp");
        written_bytes += m_comp2char.serialize(out, child, "m_comp2char");
        written_bytes += m_C.serialize(out, child, "m_C");
        written_bytes += util::write_member(m_sigma, out, child, "m_sigma");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;	
	}

    void byte_alphabet_stategy::load(std::istream& in){
		m_char2comp.load(in);
		m_comp2char.load(in);
		m_C.load(in);
		util::read_member(m_sigma, in);
	}
}
