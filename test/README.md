# Library test suite

This directory contains test code for various data structures.
We used the [googletest][GTEST] framework for testing. 

A call of `make test` in your cmake build directory will execute all tests. 
If you only want to run a test of a specific component `X` then run 
`make X`,  where X should be in the following list:

  * `bits-test` (tests basic bit operations)
  * `int-vector-test` (tests  [int_vector](../include/sdsl/int_vector.hpp))
  * `int-vector-buffer-test` (tests  [int_vector_buffer](../include/sdsl/int_vector_buffer.hpp))
  * `bit-vector-test` (tests  [bit_vector](../include/sdsl/bit_vectors.hpp) strucutres)
  * `rank-support-test` (tests  [rank_support](../include/sdsl/rank_support.hpp) structures)
  * `select-support-test` and `select-support-0-test` 
     (tests  [select_support](../include/sdsl/select_support.hpp) structures)
  * `wt-byte-test` (tests [wavelet trees](../include/sdsl/wavelet_trees.hpp) on byte alphabets)
  * `wt-int-test` (tests [wavelet trees](../include/sdsl/wavelet_trees.hpp) on integer alphabets)
  * `csa-byte-test` (tests [CSAs](../include/sdsl/suffix_arrays.hpp) on byte alphabets)
  * `csa-int-test` (tests [CSAs](../include/sdsl/suffix_arrays.hpp) on integer alphabets)
  * `cst-byte-test` (tests [CSTs](../include/sdsl/suffix_trees.hpp) on byte alphabets)
  * `cst-int-test` (tests [CSTs](../include/sdsl/suffix_trees.hpp) on integer alphabets)
  * `rmq-test` (tests [RMQ structures](../include/sdsl/rmq_support.hpp))

Test inputs are downloaded as needed before the first execution of the test.
See the [download.config](./download.config) files for details on the sources.

Executing `make test` should take about 60 minutes on a recent machine. 

Please report, if a test fails. Thanks. 

## Customization

  * Tests can be customized by editing the `.config` files.
    Test files should be located in [test_cases](./test_cases).


## Acknowledgements
  We thank 
  * Project Gutenberg for providing text files `faust.txt` and
    `zarathustra.txt`.
  * Shane Culpepper for providing the test inputs 
    `keeper.int` and `moby.int` for the integer-alphabet CSAs and CSTs.


[VG]: http://valgrind.org/ "Valgrind"
[PG]: http://www.gutenberg.org/ "Project Gutenberg"
[GTEST]: https://github.com/google/googletest "Google C++ Testing Framework"
