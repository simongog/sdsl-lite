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
/*! \file suffixtrees.hpp
    \brief suffixtrees.hpp contains generic classes for different suffix tree classes.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SUFFIXTREES
#define INCLUDED_SDSL_SUFFIXTREES

#include "sdsl_concepts.hpp"
#include "util.hpp"
#include <iostream>

using std::cout;
using std::endl;

namespace sdsl{
        
// Gets ISA[SA[idx]+d]
// d = depth of the character 0 = first position
template<class Csa>
typename Csa::size_type get_char_pos(typename Csa::size_type idx, typename Csa::size_type d, const Csa &csa){
	if( d == 0 )
		return idx;
	// if we have to apply \f$\LF\f$ or \f$\Phi\f$ more
	// than 2*d times to calc csa(csa[idx]+d), we opt to 
	// apply \f$ \Phi \f$ d times
	if( (csa.sa_sample_dens - 1) + (csa.isa_sample_dens - 1) > 2*d ){
		for(typename Csa::size_type i=0; i < d; ++i)
			idx = csa.psi[idx];
		return idx;	
	}
	return csa( csa[idx] + d );
}

// output some informations about the constructed compressed suffix tree
template<class Cst>
void cst_info(const Cst &cst){
	double csa_size, lcp_size, cst_size, nav_size, text_size;
	cout << "# alphabet size: " << (int)cst.csa.sigma << endl;
	cout << "# text size: " << (text_size = cst.csa.size()) << endl;
	cout << "# nodes: " << cst.nodes() << endl;
	cout << "# cst in MB: " << (cst_size=(((double)util::get_size_in_bytes(cst))/(1024.0*1024.0))) << endl;
	cout << "# csa in MB: " << (csa_size=(((double)util::get_size_in_bytes(cst.csa))/(1024.0*1024.0))) << endl;
	cout << "# csa samples in MB: " << ((double)util::get_size_in_bytes(cst.csa.sa_sample)+util::get_size_in_bytes(cst.csa.isa_sample))/(1024.0*1024.0) << endl;
	cout << "# lcp in MB: " << (lcp_size=(((double)util::get_size_in_bytes(cst.lcp))/(1024.0*1024.0))) << endl;
	cout << "# nav in MB: " << (nav_size=(cst_size-lcp_size-csa_size)) << endl;
	text_size /= 1024*1024;
	cout << "# "<< 8*csa_size/text_size <<" " << 8*lcp_size/text_size << " " << 8*nav_size/text_size << endl;
}

}


/** \defgroup cst Compressed Suffix Trees (CST) 
 *   This group contains data structures for compressed suffix trees. The following methods are supported:
 *    - root()
 *    - child(v,c)
 *    - ith_child(v)
 *    - ith_leaf(i)
 *    - parent(v)
 *    - sl(v)
 *    - lca(v,w)
 *    - ..
 */

#include "cst_sct.hpp"
#include "cst_sct2.hpp"
#include "cst_sct3.hpp"
#include "cst_sada.hpp"

#include "csa_uncompressed.hpp"
#include "int_vector.hpp"

#include "cst_construct.hpp"

/*
namespace sdsl{

typedef cst_sct<csa_uncompressed, int_vector<0> > cst_sct_uncompressed;
typedef cst_sada<csa_uncompressed, int_vector<0> > cst_sada_uncompressed;

}
*/

#include <iostream>
#include <string>

inline unsigned char _replace_sentinel(unsigned char c){
	if(c==0)
		return '$';
	return c;
}

inline unsigned char _vc(unsigned char c){
	if(c=='!')
		return '_';
	else
		return '!';
}

//! Output the suffix tree in tikz-format to stdout
template<class Cst>
void output_cst_in_tikz(const Cst &cst){
	std::cout<<"\
\\documentclass{article}\n\
\\usepackage{tikz}\n\
\\usepackage{verbatim}\n\
\\begin{document}\n\
\\begin{tikzpicture}\n\
[scale=0.8, transform shape, inner sep=1mm, font=\\small,\n\
innernode/.style={rectangle,draw=blue!50,fill=blue!20,thick,minimum width=#1,minimum height=0.5cm,rounded corners=2mm,anchor=south},\n\
innernode/.default=4cm,\n\
leafnode/.style={rectangle,draw=black!50,fill=black!20,thick,minimum width=#1,minimum height=0.5cm,anchor=south},\n\
leafnode/.default=1cm,\n\
]\n";

	typedef typename Cst::node_type node_type;
	typedef typename Cst::size_type size_type;
	for(typename Cst::iterator it = cst.begin(); it!=cst.end(); ++it){
		if( it.visit()==1 ){
			node_type v = *it;
			double f = 1;//1.5;
			double fy = 0.9;//1.3;
			std::string style = cst.is_leaf(v) ? "leafnode" : "innernode";
			double xpos = (cst.rb(v) + cst.lb(v));
			double ypos = cst.depth(v);
			std::cout<<"\\node["<<style<<"="<<f*cst.leaves_in_the_subtree(v)-0.2<<"cm] (node "<<ypos<<"x"<<xpos<<") at ("<<f*xpos/2<<","<<-fy*ypos<<") {"<<v<<"};"<<std::endl;
			
			if( v != cst.root() ){
				node_type p =  cst.parent(v);
				double pypos = cst.depth(p);
				std::cout<<"\\draw[->] ("<<f*xpos/2<<","<<-fy*pypos<<") -- (node "<<ypos<<"x"<<xpos<<".north);"<<std::endl;

				for(size_type i=cst.depth(p)+1; i <= cst.depth(*it); ++i){
					unsigned char c = _replace_sentinel(cst.edge(*it,i));
					double y = -fy*(pypos -1 + i-cst.depth(p) + 0.5 );
					std::cout << "\\node[anchor=south east] at ("<<f*xpos/2<<","<<y<<"){\\verb"<<_vc(c)<<c<<_vc(c)<<"};"<<std::endl;
				}
			}
		}
	}	
	std::cout<<"\
\\end{tikzpicture}\n\
\\end{document}\n";
}


#endif
