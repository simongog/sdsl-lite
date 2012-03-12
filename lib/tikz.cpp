#include "sdsl/tikz.hpp"

namespace sdsl{

void begin_tikzpicture(ostream &out, string options){ 
	out << "\\begin{tikzpicture}[" << options << "]%\n" << endl; 
}

void end_tikzpicture(ostream &out){
	out << "\\end{tikzpicture}%\n" << endl; 
}

void begin_tikzscope(ostream &out, string options){
	out << "\\begin{scope}[" << options << "]\n" << endl; 
}

void end_tikzscope(ostream &out){
   	out << "\\end{scope}%\n" << endl; 
}

void tikz_node(ostream &out, string content, string at, string name, string options){ 
	out<< "\\node[" << options << "] (" << name << ") at ("<< at <<") {" << content << "};%\n"; 
}

void tikz_coordinate(ostream &out, string at, string name, string option){ 
	out<< "\\coordinate[" << option << "] (" << name << ")";
    if( at != "" )
		out <<" at ("<< at <<")";
	out << ";%\n"; 
}

void write_y_column(ostream &out, size_t n){
	begin_tikzscope(out, "st_y");
	tikz_coordinate(out, "0,0", "y0"); 
	for(size_t i=1; i<n; ++i){
		tikz_coordinate(out, "", "y"+util::to_string(i), "below of=y"+util::to_string(i-1));
	}
	end_tikzscope(out);
}
	
}// end namespace
