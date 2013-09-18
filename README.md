SDSL - Succinct Data Structure Library
=========

What is it?
-----------

The Succinct Data Structure Library (SDSL) is a powerful and flexible C++11
library implementing succinct data structures. In total, the library contains
work of 40 [research publications](https://github.com/simongog/sdsl-lite/wiki/Literature).
Succinct data structures
can represent an object (such as a bitvector or a tree) in space close the
information-theoretic lower bound of the object while supporting operations
of the original object efficiently. The theoretical time complexity of an
operations performed on the classical data structure and the equivalent
succinct data structure are (most of the time) identical.

Why SDSL?
--------

Succinct data structures have very attractive theoretical properties. However,
in practice implementing succinct data structures is non-trivial as they are
often composed of complex operations on bitvectors. The SDSL Library provides
high quality, open source implementations of many succinct data structures
proposed in literature.

Specifically, the aim of the library is to provide basic and complex succinct
data structure which are

  * Easy and intuitive to use (like the [STL][STL], which provides classical data structures),
  * Faithful to the original theoretical results,
  * Capable of handling large inputs (yes, we support 64-bit),
  * Provide efficient construction of all implemented succinct data structures,
  while at the same time enable good run-time performance.

<img align="right" src="extras/resources/space-vis.png?raw=true" />

In addition we provide additional functionality which can help you use succinct
data structure to their full potential.

  * Each data structure can easily be serialized and loaded to/from disk.
  * We provide functionality which helps you analyze the storage requirements of any
  SDSL based data structure (see right)
  * We support features such as hugepages and tracking the memory usage of each
  SDSL data structure.
  * All implemented data structures are composable. For example, a compressed
  suffix tree can be composed of a variety different, smaller succinct data structures.
  We provide an intuitive template interface to declare the succinct data structure
  you want to use.
  * We maintain an extensive collection of examples which help you use the different
  features provided by the library.
  * All data structures are tested for correctness using a unit-testing framework.
  * We provide a large collection of supporting documentation consisting of examples,
    cheat sheets, tutorial slides and walk-through.

The library contains many succinct data structures from the following categories:

 * Bitvectors supporting Rank and Select
 * Integer Vectors
 * Wavelet Trees
 * Compressed Suffix Arrays (CSA)
 * Balanced Parentheses Representations
 * Longest Common Prefix (LCP) Arrays
 * Compressed Suffix Trees (CST)
 * Range Minimum/Maximum Query (RMQ) Structures

For a complete overview including theoretical bounds see the
[cheat sheet](extras/cheatsheet/sdsl-cheatsheet.pdf?raw=true) or the
[wiki](https://github.com/simongog/sdsl-lite/wiki/List-of-Implemented-Data-Structures).

Documentation
-------------

We provide an extensive set of documentation describing all data structures
and features provided by the library. Specifically we provide

* A [cheat sheet](extras/cheatsheet/sdsl-cheatsheet.pdf?raw=true) which succinctly
describes the usage of the library.
* A set of [example](examples/) programs demonstrating how different features
of the library are used.
* A tutorial [presentation][TUT] with the [example code](tutorial/) using in the
sides demonstrating all features of the library in a step-by-step walk-through.
* [Unit Tests](test/) which contain small code snippets used to test each
library feature.

Requirements
------------

The SDSL library requires:

* A modern, C++11 ready compiler such as `g++` version 4.7 or higher or `clang` version 3.2 or higher.
* The [cmake][cmake] build system.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* For increased performance the processor of the system should support fast bit operations available in `SSE4.2`

Installation
------------

To download and install the library use the following commands.

```sh
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
```

This installs the sdsl library into the `include` and `lib` directories in your
home directory. A different location prefix can be specified as a parameter of
the `install.sh` script:

```sh
./install /usr/local/
```

To remove the library from your system use the provided uninstall script:

```sh
./uninstall.sh
```

Getting Started
------------

To get you started with the library you can start by compiling the following
sample program which constructs a compressed suffix array (a FM-Index) over the
text `mississippi!`, counts the number of occurrences of pattern `si` and
stores the data structure, and a space usage visualization to the
files `fm_index-file.sdsl` and `fm_index-file.sdsl.html`:

```cpp
#include <sdsl/suffix_arrays.hpp>

int main() {
  sdsl::csa_wt<> fm_index;
  sdsl::construct_im(fm_index, "mississippi!", 1);
  std::cout << "'si' occurs " << sdsl::count(fm_index,"si") << " times.\n";
  sdsl::store_to_file(fm_index,"fm_index-file.sdsl");
  sdsl::write_structure<HTML_FORMAT>(fm_index,"fm_index-file.sdsl.html");
}
```

To compile the program using `g++` run:

```sh
g++ -std=c++11 -O3 -I ~/include -L ~/lib program.cpp -o lsdsl
```

Next we suggest you look at the comprehensive [tutorial][TUT] of Simon Gog which describes
all major features of the library or look at some of the provided [examples](examples).

Test
----

Implementing succinct data structures can be tricky. To ensure that all data
structures behave as expected, we created a large collection of unit tests
which can be used to check the correctness of the library on your computer.
The [test](./test) directory contains test code. We use [googletest][GTEST]
framework and [make][MAKE] to run the tests. See the README file in the
directory for details.

To simply run all unit tests type

```sh
cd sdsl-lite/test
make
```

Note: Running the tests requires several sample files to be downloaded from the web
and can take up to 2 hours on slow machines.


Benchmarks
----------

To ensure the library runs efficiently on your system we suggest you run our
[benchmark suite](benchmark). The benchmark suite recreates a
popular [experimental study](http://arxiv.org/abs/0712.3360) which you can
directly compare to the results of your benchmark run.

Bug Reporting
------------

While we use an extensive set of unit tests and test coverage tools you might
still find bugs in the library. We encourage you to report any problems with
the library via the [github issue tracking system](https://github.com/simongog/sdsl-lite/issues)
of the project.

The Latest Version
------------------

The latest version can be found on the SDSL github project page https://github.com/simongog/sdsl-lite .

If you are running experiments in an academic settings we suggest you use the
most recent [released](https://github.com/mpetri/sdsl-lite/releases) version
of the library. This allows others to reproduce your experiments exactly.

Licensing
---------

The SDSL library is provided under the GNU General Public License (GPLv3). For
more information see the COPYING file in the library directory.

Lots of time was spent implementing the many features of the library. If you
use the library in an academic setting please cite the following paper:

_Simon Gog, Matthias Petri: Optimized Succinct Data Structures for Massive Data, Accepted for publication in Software, Practice and Experience_.

## External Resources used in SDSL

We have included the code of two excellent suffix array
construction algorithms.

* Yuta Mori's incredible fast suffix [libdivsufsort][DIVSUF]
  algorithm (version 2.0.1) for byte-alphabets.
* An adapted version of Jesper Larsson's implementation of the
  algorithm of [Larson and Sadakane][LS] for integer-alphabets.

Additionally, we use the [googletest][GTEST] framework to provide unit tests.
Our visualizations are implemented using the [d3js][d3js]-library.

Authors
--------

The main contributors to the library are:

* [Timo Beller](https://github.com/tb38)
* [Simon Gog](https://github.com/simongog) (Creator)
* [Matthias Petri](https://github.com/mpetri)

This project further profited from excellent input of many coders. Stefan
Arnold helped us with tricky template questions. We are also grateful to Kalle Karhu,
Dominik Kempa, and Shanika Kuruppu for bug reports.

Contribute
----------

Are you working on a new or improved implementation of a succinct data structure?
We encourage you to contribute your implementation to the SDSL library to make
your work accessible to the community within the existing library framework.
Feel free to contact any of the authors or create an issue on the
[issue tracking system](https://github.com/simongog/sdsl-lite/issues).


[STL]: http://www.sgi.com/tech/stl/ "Standard Template Library"
[pz]: http://pizzachili.di.unipi.it/ "Pizza&amp;Chli"
[d3js]: http://d3js.org "D3JS library"
[cmake]: http://www.cmake.org/ "CMake tool"
[MAKE]: http://www.gnu.org/software/make/ "GNU Make"
[gcc]: http://gcc.gnu.org/ "GNU Compiler Collection"
[SPE]: http://people.eng.unimelb.edu.au/sgog/optimized.pdf "Preprint SP&amp;E article"
[DIVSUF]: http://code.google.com/p/libdivsufsort/ "libdivsufsort"
[LS]: http://www.sciencedirect.com/science/article/pii/S0304397507005257 "Larson &amp; Sadakane Algorithm"
[GTEST]: https://code.google.com/p/googletest/ "Google C++ Testing Framework"
