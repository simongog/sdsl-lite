SDSL - Succinct Data Structure Library
=========

[![Build Status](https://travis-ci.org/simongog/sdsl-lite.svg?branch=master)](https://travis-ci.org/simongog/sdsl-lite)

What is it?
-----------

The Succinct Data Structure Library (SDSL) is a powerful and flexible C++11
library implementing succinct data structures. In total, the library contains
the highlights of 40 [research publications][SDSLLIT]. Succinct data structures
can represent an object (such as a bitvector or a tree) in space close to the
information-theoretic lower bound of the object while supporting operations
of the original object efficiently. The theoretical time complexity of an
operation performed on the classical data structure and the equivalent
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

<a href="http://simongog.github.io/assets/data/space-vis.html" >
<img align="right" src="extras/resources/space-vis.png?raw=true" />
</a>

In addition we provide additional functionality which can help you use succinct
data structure to their full potential.

  * Each data structure can easily be serialized and loaded to/from disk.
  * We provide functionality which helps you analyze the storage requirements of any
  SDSL based data structure (see right)
  * We support features such as hugepages and tracking the memory usage of each
  SDSL data structure.
  * Complex structures can be configured by template parameters and therefore
    easily be composed. There exists one simple method which constructs
    all complex structures.
  * We maintain an extensive collection of examples which help you use the different
    features provided by the library.
  * All data structures are tested for correctness using a unit-testing framework.
  * We provide a large collection of supporting documentation consisting of examples,
    [cheat sheet][SDSLCS], [tutorial slides and walk-through][TUT].

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
[cheat sheet][SDSLCS] or the
[wiki](https://github.com/simongog/sdsl-lite/wiki/List-of-Implemented-Data-Structures).

Documentation
-------------

We provide an extensive set of documentation describing all data structures
and features provided by the library. Specifically we provide

* A [cheat sheet][SDSLCS] which succinctly
describes the usage of the library.
* An doxygen generated [API reference][DOXYGENDOCS] which lists all types and functions
of the library.
* A set of [example](examples/) programs demonstrating how different features
of the library are used.
* A tutorial [presentation][TUT] with the [example code](tutorial/) using in the
sides demonstrating all features of the library in a step-by-step walk-through.
* [Unit Tests](test/) which contain small code snippets used to test each
library feature.

Requirements
------------

The SDSL library requires:

* A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
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
./install.sh /usr/local/
```

To build a portable sdsl library without using `SSE4.2` and `AVX2` instructions, set `BUILD_PORTABLE` at build time, e.g. `BUILD_PORTABLE=1 ./install.sh` or `mkdir build && BUILD_PORTABLE=1 cmake ..`.
These instructions are enabled by default if the processor of the build system supports them.

To remove the library from your system use the provided uninstall script:

```sh
./uninstall.sh
```

There is also a [Gentoo Ebuild for SDSL][EBUILD] by Mathias Weller.

Getting Started
------------

To get you started with the library you can start by compiling the following
sample program which constructs a compressed suffix array (a FM-Index) over the
text `mississippi!`, counts the number of occurrences of pattern `si` and
stores the data structure, and a space usage visualization to the
files `fm_index-file.sdsl` and `fm_index-file.sdsl.html`:

```cpp
#include <sdsl/suffix_arrays.hpp>
#include <fstream>

using namespace sdsl;

int main() {
    csa_wt<> fm_index;
    construct_im(fm_index, "mississippi!", 1);
    std::cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
    store_to_file(fm_index,"fm_index-file.sdsl");
    std::ofstream out("fm_index-file.sdsl.html");
    write_structure<HTML_FORMAT>(fm_index,out);
}
```

To compile the program using `g++` run:

```sh
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib program.cpp -o program -lsdsl -ldivsufsort -ldivsufsort64
```

Next we suggest you look at the comprehensive [tutorial][TUT] which describes
all major features of the library or look at some of the provided [examples](examples).

Test
----

Implementing succinct data structures can be tricky. To ensure that all data
structures behave as expected, we created a large collection of unit tests
which can be used to check the correctness of the library on your computer.
The [test](./test) directory contains test code. We use [googletest][GTEST]
framework and [make][MAKE] to run the tests. See the README file in the
directory for details.

To simply run all unit tests after installing the library type

```sh
cd sdsl-lite/build
make test-sdsl
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
most recent [released](https://github.com/simongog/sdsl-lite/releases) version
of the library. This allows others to reproduce your experiments exactly.

Licensing
---------

The SDSL library is free software provided under the GNU General Public License
(GPLv3). For more information see the [COPYING file][CF] in the library
directory.

We distribute this library freely to foster the use and development of advanced
data structure. If you use the library in an academic setting please cite the
following paper:

    @inproceedings{gbmp2014sea,
      title     = {From Theory to Practice: Plug and Play with Succinct Data Structures},
      author    = {Gog, Simon and Beller, Timo and Moffat, Alistair and Petri, Matthias},
      booktitle = {13th International Symposium on Experimental Algorithms, (SEA 2014)},
      year      = {2014},
      pages     = {326-337},
      ee        = {http://dx.doi.org/10.1007/978-3-319-07959-2_28}
    }

A preliminary version is available [here on arxiv][SEAPAPER].

## External Resources used in SDSL

We have included the code of two excellent suffix array
construction algorithms.

* Yuta Mori's incredible fast suffix [libdivsufsort][DIVSUF]
  algorithm for byte-alphabets.
* An adapted version of [Jesper Larsson's][JESL] [implementation][QSUFIMPL] of
  suffix array sorting on integer-alphabets (description of [Larsson and Sadakane][LS]).

Additionally, we use the [googletest][GTEST] framework to provide unit tests.
Our visualizations are implemented using the [d3js][d3js]-library.


Authors
--------

The main contributors to the library are:

* [Johannes Bader](https://github.com/olydis)
* [Timo Beller](https://github.com/tb38)
* [Simon Gog](https://github.com/simongog) (Creator)
* [Matthias Petri](https://github.com/mpetri)

This project is also supported by code contributions
from other researchers. E.g. Juha Kärkkäinen,
[Dominik Kempa](https://github.com/dkempa),
and Simon Puglisi contributed a compressed bitvector
implementation ([hyb_vector][HB]).
This project further profited from excellent input of our students
Markus Brenner, Alexander Diehm, Christian Ocker, and Maike Zwerger. Stefan
Arnold helped us with tricky template questions. We are also grateful to
[Diego Caro](https://github.com/diegocaro),
[Travis Gagie](https://github.com/TravisGagie),
Kalle Karhu,
[Bruce Kuo](https://github.com/bruce3557),
Jan Kurrus,
[Shanika Kuruppu](https://github.com/skuruppu),
Jouni Siren,
and [Julio Vizcaino](https://github.com/garviz)
for bug reports.

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
[DIVSUF]: https://github.com/y-256/libdivsufsort/ "libdivsufsort"
[LS]: http://www.sciencedirect.com/science/article/pii/S0304397507005257 "Larson &amp; Sadakane Algorithm"
[GTEST]: https://code.google.com/p/googletest/ "Google C++ Testing Framework"
[SDSLCS]: http://simongog.github.io/assets/data/sdsl-cheatsheet.pdf "SDSL Cheat Sheet"
[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature "Succinct Data Structure Literature"
[TUT]: http://simongog.github.io/assets/data/sdsl-slides/tutorial "Tutorial"
[QSUFIMPL]: http://www.larsson.dogma.net/qsufsort.c "Original Qsufsort Implementation"
[JESL]: http://www.itu.dk/people/jesl/ "Homepage of Jesper Larsson"
[CF]: https://github.com/simongog/sdsl-lite/blob/master/COPYING "Licence"
[SEAPAPER]: http://arxiv.org/pdf/1311.1249v1.pdf "SDSL paper"
[HB]: https://github.com/simongog/sdsl-lite/blob/hybrid_bitvector/include/sdsl/hybrid_vector.hpp "Hybrid bitevctor"
[DOXYGENDOCS]: http://algo2.iti.kit.edu/gog/docs/html/index.html "API Reference"
[EBUILD]: https://github.com/igel-kun/silly_ebuilds "Gentoo EBuild"
