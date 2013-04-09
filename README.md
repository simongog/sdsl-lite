# SDSL: Succinct Data Structure Library

This is a C++ template library for succinct data structures
called _sdsl_. 

Succinct data structures are fascinating: They represent an
object (like a bitvector, a tree, suffix array,...) in space
close the information-theoretic lower bound of the object
but the defined operations can still be performed 
efficiently. Hmmm, at least in theory ;) Actually there
is still a big gap between theory and practice. Why?
The time complexity of an operations performed
on the classical fat data structure and the slim 
succinct data structure are the same most time in
theory. However, in practice succinct structures are
slow since the operations require often memory accesses
with bad locality of references. Moreover, often the in theory
small sub-linear space data structures account for
a large amount of memory, since they are only 
asymptotic sub-linear and the input size for which
they are negligible in practice is galactic.

The aim of the library is to provide basic and complex succinct
data structure which are
  * easy to use (like the [STL][STL], which provides classical data structures),
  * capable of handling large inputs (yes, we support 64-bit),
  * provide excellent performance in construction, and 
  * provide excellent operation performance 

We used several techniques to reach the performance goals:
  * bit-parallelism on 64-bit word
  * cache-friendly data structure layouts
  * semi-external construction algorithms
  * the built-in popcount operation
  * hugepages to bypass the TLB bottleneck

Details are in our [comprehensive experimental study][SPE].

## List of implemented data structures

 * Bitvectors
   * An uncompressed mutual bitvector (`bit_vector`)
   * An uncompressed immutable bitvector ([bit_vector_il](./include/sdsl/bit_vector_il.hpp))
   * A ![H_0](http://latex.codecogs.com/gif.latex?H_0)-compressed
     immutable bitvector ([rrr_vector<>](./include/sdsl/rrr_vector.hpp))
   * A bitvector for sparse populated arrays ([sd_vector<>](./include/sdsl/sd_vector.hpp))
 * Rank Support (RS) and Select Support (SS)
   * Several rank and select implementations with different time-space
     trade-offs for the uncompressed bitvectors 
     ([rank_support_v](./include/sdsl/rank_support_v.hpp),
      [rank_support_v5](./include/sdsl/ran_support_v5.hpp),
      [select_support_mcl](./include/sdsl/select_support_mcl.hpp), ...)
   * Rank and select for compressed bitvectors 
     ([rank_support_rrr](./include/sdsl/rrr_vector.hpp),
      [rank_support_sd](./include/sdsl/sd_vector.hpp), ...)
 * Variable-length Coders
   * Elias- ![\delta](http://latex.codecogs.com/gif.latex?%5Cdelta) coder 
     ([coder::elias_delta](./include/sdsl/coder_elias_delta.hpp))
   * Fibonacci-coder ([coder::fibonacci](./include/sdsl/coder_fibonacci.hpp))
 * Integer Vectors
   * Mutable vectors for (compile-time) fixed `w`-bit integers ([int_vector<w>](./include/sdsl/int_vector.hpp))
   * Mutable vector for (run-time) fixed `w`-bit integers ([int_vector<0>](./include/sdsl/int_vector.hpp),
     `w` passed to the constructor)
   * Immutable compressed integer vectors using a variable-length coder `coder` 
     ([enc_vector<coder>](./include/sdsl/enc_vector.hpp), 
     [vlc_vector<coder>](./include/sdsl/vlc_vector.hpp))
 * Wavelet Trees (WT) (all immutable)
   * Balanced wavelet tree for a byte-alphabet ([wt](./include/sdsl/wt.hpp))
   * Balanced wavelet tree for a integer-alphabet ([wt_int](./include/sdsl/wt_int.hpp))
   * Huffman-shaped wavelet tree for a byte-alphabet ([wt_huff](./include/sdsl/wt_huff.hpp)) 
   * Run-length compressed wavelet trees for a byte-alphabet 
     ([wt_rlmn](./include/sdsl/wt_rlmn.hpp),
      [wt_rlg](./include/sdsl/wt_rlg.hpp),
      [wt_rlg8](./include/sdsl/wt_rlg8.hpp))
 * Compressed Suffix Arrays (CSA) (all immutable)
   * [csa_bitcompressed](./include/sdsl/csa_bitcompressed.hpp) is based on the bitcompressed SA and inverse SA.
   * [csa_wt](./include/sdsl/csa_wt.hpp) is based on a WT of the BWT.
   * [csa_sada](./include/sdsl/csa_sada.hpp) is based on the compressed
      ![\Psi](http://latex.codecogs.com/gif.latex?%5CPsi)-function
 * Balanced Parentheses Support (BPS) (all immutable)
   * A range-min-max-tree implementation ([bp_support_sada](./include/sdsl/bp_support_sada.hpp))
     to support operations `find_open`, `find_close`, `enclose`, `double_enclose`,...
   * Hierarchical solutions with pioneer parentheses 
     ([bp_support_g](./include/sdsl/bp_support_g.hpp),
      [bp_support_gg](./include/sdsl/bp_support_gg.hpp))
 * Longest Common Prefix (LCP) Arrays (all immutable)
   * [lcp_bitcompressed](./include/sdsl/lcp_bitcompressed.hpp) is a bitcompressed version
   * [lcp_byte](./include/sdsl/lcp_byte.hpp) encodes small values with one
     byte and large values with two words 
   * [lcp_dac](./include/sdsl/lcp_dac.hpp) used direct accessible codes 
   * [lcp_wt](./include/sdsl/lcp_wt.hpp) stores small values in a WT and
     large value in on word.
   * [lcp_vlc](./include/sdsl/lcp_vlc.hpp) stores the values in a `vlc_vector`.
   * [lcp_support_sada](./include/sdsl/lcp_support_sada.hpp) uses a bitvector of 2n bits, a
     select structure supporting it, and the corresponding CSA.
   * [lcp_support_tree](./include/sdsl/lcp_support_tree.hpp) uses the topology of the
     corresponding CST.
   * [lcp_support_tree2](./include/sdsl/lcp_support_tree2.hpp) uses the corresponding CSA 
	 and CST.
 * Compressed Suffix Trees (CST) (all immutable)
   * [cst_sada](./include/sdsl/cst_sada.hpp) provides very fast navigation operations; worst case  
     space |CSA|+|LCP|+4n+o(n)
   * [cst_sct3](./include/sdsl/cst_sct3.hpp) representing nodes as intervals in the suffix array;
     worst case space |CSA|+|LCP|+3n+o(n)
 * Range Minimum/Maximum Query (RMQ) Structures (all immutable)
   * Self-contained RMQ structure using 2n+o(n) bits or 4n+o(n) bits
     ([rmq_succinct_sct](./include/sdsl/rmq_succinct_sct.hpp),
      [rmq_succinct_sada](./include/sdsl/rmq_succinct_sada.hpp))
   * Non-succinct support structure for RMQ ([rmq_support_sparse_table](./include/sdsl/rmq_support_sparse_table.hpp))


## Example of a complex data structure

Let us now show how you can assemble even a very
complex data structure very easily. Lets begin with
the most complex one, a CST!
It basically consists of a CSA, an compressed LCP-array,
and a succinct representation of the tree topology;
each part can be specified by a template parameter.
Say, we want fast navigation operations, so we take
the class `cst_sada<cst_type, lcp_type, bp_support_type>` 
for our CST. Now we can specify the type of CSA.
Lets take a CSA based on wavelet tree:
`csa_wt<wt_type, SA_sample_dens, inv_SA_sample_dens>`.
We can recursively specify the used types. So
now we can specify the used wavelet tree, say
a run-length compressed wavelet tree 
(`wt_rlmn<>`). We could recurse again and specify, each detail
of the wavelet tree (e.g. which rank support structure 
should be used) but we stick now with the default
configuration which uses an `sd_vector` for the
marking of the heads of the runs in the wavelet tree.
Lets choose at last a LCP array which uses
the topology of the CST and the CSA to
compress the LCP values (`lcp_support_tree2`) and
stick with default template parameters for all
types. So the final type looks like this:
``cst_sada<cst_wt<wt_rlmn<> >, lcp_support_tree2<> >``.

Now, lets explore the data structure a little bit.
We take the english.100MB input from 
the [Pizza&amp;Chili][pz]-corpus,
construct the CST-object, output its structure, and visualise
it using the [d3js][d3js]-library. Have fun with 
[the result][DBLPCSTRES]


## Types of data structures

The data structures in the library can be divided into several classes:
  * Objects of _mutable_ classes can be changed after construction (e.g.
    we can assign new values to the elements of an `int_vector`)
  * Objects of _immutable_ classes can not be changed after construction 
    (e.g. you can not assign a new value to an element of a
    compressed suffix array, say `csa_wt`)
  * Objects of _support_ classes add functionality to objects of
    self-contained classes. For example an object of type `rank_support_v`
    addes constant time `rank(i)`-functionality to an object of type
    `bit_vector`, or an object of of type `bp_support_sada` adds
    `find_open(i)`-functionality to a `bit_vector` object, which
    represents a balanced parentheses sequence.

Each _sdsl_-class `X` has to implement the following methods:
  * The standard constructor `X()`
  * The copy constructor `X(const &X)`
  * Swap operator `swap(const &X)`
  * serialize operator `serialize(std::ostream &out, structure_tree_node* v, std::string name)`
  * load operator `load(std::istream &in)`

We provide many handy methods for _sdsl_ objects: 
  * `store_to_file(const T &o, std::string file)` stores object `o` to `file`
  * `load_from_file(T &o, std::string file)` loads object `o` from `file`
  * `size_in_bytes(const T &o)` returns the number of bytes needed to represent 
     object `o` in memory.
  * `write_structure<FORMAT>(const T &o, std::ostream &out)` writes the structure
     of object `o` in JSON or R format (`FORMAT`=`JSON_FORMAT` or `R_FORMAT`) into `out`
  * `util::clear(T &o)` frees space by setting o to the empty object.
  * ...for more have a look into the cheat sheet in `extras/cheatsheet`.

## Supported platforms

The library was successfully tested on the following
configurations
  * Mac OS X 10.7.3 on a MacBookPro equipped with a Intel Core i5 CPU
  * Ubuntu Linux 12.04 running on a server equipped with Intel Xeon (E5640) CPUs


## Installation

The installation requires that the [cmake tool][cmake]
and a C++ compiler (e.g. from the [GNU Compiler Collection][gcc]
is installed.
You can than install the library into an directory `SDSL_INSTALL_DIR` by
calling
``
./install SDSL_INSTALL_DIR
``
If `SDSL_INSTALL_DIR` is not specified your home directory is used.
The library header files will be located in the directory
`SDSL_INSTALL_DIR/include` and the library in the directory
`SDSL_INSTALL_DIR/lib`. After the installation you can 
execute the tests in the `test` directory or start
with some code examples in the `examples` folder.
A cheat sheet can be generated my running make in the 
`extras/cheatsheet` folder.


## Tests

The [test](./test) directory contains test code for many library
structures. We use [googletest][GTEST] framework and [make][MAKE]
to run the tests. See the README file in the directory for details.

## Benchmarks

You can find out how efficient the library works on your system by
running experiments in the [benchmark](./benchmark) directory.

## Examples

Compile the examples with `make` and experience
how easy it is to use succinct data structures.


## Construction of Suffix Arrays

We have included the code of two excellent suffix array
construction algorithms.

* Yuta Mori's incredible fast suffix [libdivsufsort][DIVSUF]
  algorithm (version 2.0.1) for byte-alphabets.
* An adapted version of Jesper Larsson's implementation of the
  algorithm of [Larson and Sadakane][LS] for integer-alphabets.


## Contributors

This project profits from excellent input of many coders. Timo Beller improved the construction
process during the last month and Matthias Petri contributed new bitvectors and helped a lot
in making the library more accessible. Stefan Arnold helped us with tricky template questions.
We a grateful to Kalle Karhu, Dominik Kempa, and Shanika Kuruppu for bug reports.

## References

The header file of each data structure should contain all relevant references for the
data structure. If you find a header file, where a reference is missing, please contact
Simon.


The latest paper about the library itself is:


 * Simon Gog, Matthias Petri: _Optimized Succinct Data Structures for Massive Data_,
   Accepted for publication in Software, Practice and Experience.


[STL]: http://www.sgi.com/tech/stl/ "Standart Template Library"
[pz]: http://pizzachili.di.unipi.it/ "Pizza&amp;Chli"
[d3js]: http://d3js.org "D3JS library"
[cmake]: http://www.cmake.org/ "CMake tool"
[MAKE]: http://www.gnu.org/software/make/ "GNU Make"
[gcc]: http://gcc.gnu.org/ "GNU Compiler Collection"
[DBLPCSTRES]: http://people.eng.unimelb.edu.au/sgog/sdsl_explore/dblp.xml.100MB_cst_sada_wt_rlmn_lcp_tree2.html "CST visualization"
[SPE]: http://people.eng.unimelb.edu.au/sgog/optimized.pdf "Preprint SP&amp;E article"
[DIVSUF]: http://code.google.com/p/libdivsufsort/ "libdivsufsort"
[LS]: http://www.sciencedirect.com/science/article/pii/S0304397507005257 "Larson &amp; Sadakane Algorithm"
[GTEST]: https://code.google.com/p/googletest/ "Google C++ Testing Framework"
