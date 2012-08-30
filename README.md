SDSL: Succinct Data Structure Library
=====================================

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
  * easy to use (the library is structure like the 
    [STL](http://www.sgi.com/tech/stl/), which provides
    classical data structures),
  * capable of handling large inputs (yes, we support 64-bit),
  * provide excellent performance in construction, and 
  * provide excellent operation performance 

A lot of engineering tricks had to be applied to
reach the performance goal, for instance the use a semi-external
algorithm, bit-parallelism on 64-bit words,
and cache-friendly algorithms. 

List of implemented data structures
-----------------------------------
  * Bitvectors
    * An uncompressed mutual bitvector (`bit_vector`)
    * An uncompressed immutable bitvector (`bit_vector_interleaved`)
    * A ![H_0](http://latex.codecogs.com/gif.latex?H_0)-compressed immutable bitvector (`rrr_vector<>`)
    * A bitvector for sparse populated arrays (`sd_vector<>`)
  * Rank and Select Support Structures
    * Several rank and select implementations with different time-space
      trade-offs for the uncompressed bitvectors 
      (`rank_support_v`,`rank_support_v5`,`select_support_mcl`,...)
    * Rank and select for compressed bitvectors (`rrr_rank_support<>`, `sd_rank_support<>`,...) 
  * Variable-length Coders
    * Elias- ![\delta](http://latex.codecogs.com/gif.latex?%5Cdelta) coder (`coder::elias_delta`)
    * Fibonacci-coder (`coder::fibonacci`)
  * Integer Vectors
    * Mutable vectors for (compile-time) fixed `w`-bit integers (`int_vector<w>`)
    * Mutable vector for (run-time) fixed `w`-bit integers (`int_vector<0>`, `w` passed to the constructor)
    * Immutable compressed integer vector using a variable-length coder `coder` (`enc_vector<coder>`)
  * Wavelet Trees (all immutable)
    * Balanced wavelet tree for a byte-alphabet (`wt`)
    * Balanced wavelet tree for a integer-alphabet (`wt_int`)
    * Huffman-shaped wavelet tree for a byte-alphabet (`wt_huff`) 
    * Run-length compressed wavelet trees for a byte-alphabet (`wt_rlmn`, `wt_rlg`, and `wt_rlg8`)
  * Compressed Suffix Arrays (CSA) (all immutable)
    * CSA based on a wavelet tree (`csa_wt`)
    * CSA based on the compressed ![\Psi](http://latex.codecogs.com/gif.latex?%5CPsi)-function `csa_sada`
  * Balanced Parentheses Support Structures (all immutable)
    * A range-min-max-tree implementation (`bp_support_sada`) to `find_open`, `find_close`,
	  `enclose`, `double_enclose`,...
    * Hierarchical solution with pioneer parentheses (`bp_support_g`, `bp_support_gg`)
  * Range Minimum Support (RMQ) Structures (all immutable)
    * Self-contained RMQ structure using 2n+o(n) bits or 4n+o(n) bits (`rmq_succinct_sct`, `rmq_succinct_sada`)
    * Non-succinct support structure for RMQ (`rmq_support_sparse_table`)
  * Longest Common Prefix (LCP) Arrays (all immutable)
    * LCP-array based on direct accessible codes (`lcp_dac`)
    * LCP-array encodes small values with a byte and large values with a word (`lcp_kurtz`)
    * LCP-array encodes all values in a wavelet tree (`lcp_wt`)
    * Compressed LCP-array dependent on the corresponding CSA (`lcp_support_sada`)
    * Compressed LCP-array dependent on the corresponding CST (`lcp_support_tree`)
    * Compressed LCP-array dependent on the corresponding CSA and CST (`lcp_support_tree2`)
  * Compressed Suffix Trees(CSTs) (all immutable)
    * CST providing very fast navigation operations (`cst_sada`)
    * CST representing nodes as intervals in the suffix array (`cst_sct3`)

Example of a complex data structure
-----------------------------------
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
the [Pizza&amp;Chili](http://pizzachili.di.unipi.it/)-corpus,
construct the CST-object, output its structure, and visualise
it using the [d3js](http://d3js.org/)-library. Have fun with 
[the result](http://people.eng.unimelb.edu.au/sgog/sdsl_explore/dblp.xml.100MB_cst_sada_wt_rlmn_lcp_tree2.html).


Types of data structures
------------------------
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

We provide many handy methods for _sdsl_ objects in the `util` namespace: 
  * `util::store_to_file(const X &x, const char* file_name)` stores the object `x` to the file
  * `util::clear(X &x)` deletes the object and frees the space 
  * `util::load_from_file(X &x, const char* file_name)` loads the object `x` from the file
  * `util::assign(X &x, Y &y)` if the type of `X` equals `Y`, then `x` and `y` are swapped,
     otherwise `y` is assigned to `x` by `x = T(y)`
  * `util::get_size_in_bytes(const X &x)` returns the number of bytes needed to represent 
     object `x` in memory.
  * `util::write_structure<FORMAT>(const X &x, std::ostream &out)` writes the structure
     of the data structure in JSON or R format (`FORMAT`=`JSON_FORMAT` or `R_FORMAT`)

Supported platforms
-------------------
The library was successfully tested on the following
configurations
  * Mac OS X 10.7.3 on a MacBookPro equipped with a Intel Core i5 CPU
  * Ubuntu Linux 12.04 running on a server equipped with Intel Xeon (E5640) CPUs

We plan to support Windows in the near future.

Installation
------------
The installation requires that the [cmake tool](http://www.cmake.org/)
and a C++ compiler (e.g. from the [GNU Compiler Collection](http://gcc.gnu.org/))
is installed.
You can than install the library into a directory `SDSL_INSTALL_DIR` by
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

Tests
-----
We have used the gtest framework for the tests.
Compile with `make` and run tests with `make test`. We have another
target `vtest` which runs the test with the valgrind tool.
`make test` will try to download some texts from a
gutenberg.org mirror. See the README file in the directory for details.

Examples
--------
Compile the examples with `make` and experience
how esay it is to use succinct data structures.

Construction of Suffix Arrays
-----------------------------
The current version includes Yuta Mori's incredible fast suffix array
construction library [libdivsufsort](http://code.google.com/p/libdivsufsort/)
version 2.0.1.



Contributors
------------

Here is a list of contributes:

Code:
  * Stefan Arnold
  * Timo Beller
  * Simon Gog
  * Shanika Kuruppu
  * Matthias Petri

Bug reports:
  * Kalle Karhu
  * Dominik Kempa

New contributors are welcome any time!

Have fun with the library!

