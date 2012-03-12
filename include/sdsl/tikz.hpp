//! Author: Simon Gog (simon.gog@unimelb.edu.au)

#ifndef INCLUDED_SDSL_TIKZ
#define INCLUDED_SDSL_TIKZ

#include <string>
#include <ostream>
#include <iostream>
#include <sdsl/util.hpp>

using std::string;
using std::ostream;
using std::endl;

namespace sdsl{

void begin_tikzpicture(ostream &out, string options="");
void end_tikzpicture(ostream &out);
void begin_tikzscope(ostream &out, string options="");
void end_tikzscope(ostream &out);
void tikz_node(ostream &out, string content="", string at="0,0", string name="", string options=""); 
void tikz_coordinate(ostream &out, string at="0,0", string name="", string option="");

template<class T>
void write_tikz_column_from_container(ostream &out, const T& vec, string name_prefix="i"){
	tikz_node(out, "", "0,0cm", name_prefix, "st_"+name_prefix);
	for(typename T::size_type i=0; i < vec.size(); ++i){
		tikz_node(out, util::to_latex_string(vec[i]), name_prefix + "|- y" + util::to_string(i) , 
				  name_prefix+util::to_string(i), "st_elem_"+name_prefix);
	}
}

template<class tContainer>
void write_tikz_array(ostream &out, const tContainer &v, string array_name="", bool escape=false){
	if ( array_name != "" ){
		out << "\\def\\" << array_name << "{%" << endl;
	}else{
		out << "{";
	}
	for (typename tContainer::size_type i=0; i < v.size(); ++i) {
		if ( i > 0 )
			out << ",";
		string w = util::to_latex_string( v[i] );
		if ( escape ){
			if( w.size() > 0  and w[0]=='\\' )
				w = "\\noexpand"+w;
			out << "\"" << w << "\"";
		}else{
			out << w;
		}
	}
	out << "}";
   	if ( array_name != "" ){
		out << "%" << endl;
	}
}

void write_y_column(ostream &out, size_t n);

} // end namespace
#endif
