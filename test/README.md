# Library test suite

This directory contains test code for various data structures.
We used the [googletest][GTEST] framework for testing.

A call of `make test-sdsl` in your cmake build directory will execute all tests.
If you only want to run a test of a specific component `X` then run
`make X`,  where X should be in the following list:

  * `bits-test` (tests basic bit operations)
  * `bit-vector-test` (tests [bit_vector](../include/sdsl/bit_vectors.hpp) structures)
  * `coder-test` (tests [coder](../include/sdsl/coder.hpp) e.g. elias-gamma, fibonacci, comma)
  * `compile-test` (tests if a program that include all header-files compiles)
  * `csa-byte-test` (tests [CSAs](../include/sdsl/suffix_arrays.hpp) on byte alphabets)
  * `csa-int-test` (tests [CSAs](../include/sdsl/suffix_arrays.hpp) on integer alphabets)
  * `cst-byte-test` (tests [CSTs](../include/sdsl/suffix_trees.hpp) on byte alphabets)
  * `cst-int-test` (tests [CSTs](../include/sdsl/suffix_trees.hpp) on integer alphabets)
  * `int-vector-buffer-test` (tests [int_vector_buffer](../include/sdsl/int_vector_buffer.hpp))
  * `int-vector-mapper-test` (tests [int_vector_mapper](../include/sdsl/int_vector_mapper.hpp))
  * `int-vector-test` (tests [int_vector](../include/sdsl/int_vector.hpp))
  * `inv-perm-support-test` (tests [inv_perm_support](../include/sdsl/inv_perm_support.hpp))
  * `k2-treap-test` (tests [k2-treap](../include/sdsl/k2_treap.hpp))
  * `lcp-construct-test` (tests different lcp-array construction algorithms)
  * `nn-dict-dynamic-test` (tests [nn-dict-dynamic](../include/sdsl/nn_dict_dynamic.hpp))
  * `rank-support-test` (tests  [rank_support](../include/sdsl/rank_support.hpp) structures)
  * `rmq-test` (tests [RMQ structures](../include/sdsl/rmq_support.hpp))
  * `sa-construct-test` (tests different suffix-array construction algorithms)
  * `sd-vector-test` (tests [sd_vector](../include/sdsl/sd_vector.hpp) constructors)
  * `search-bidirectional-test` (tests bidirectional search)
  * `select-support-test` and `select-support-0-test`
     (tests [select_support](../include/sdsl/select_support.hpp) structures)
  * `sorted-int-stack-test` (tests [sorted-int-stack](../include/sdsl/sorted_int_stack.hpp))
  * `sorted-stack-support-test` (tests [sorted-stack-support](../include/sdsl/sorted_stack_support.hpp))
  * `wt-byte-test` (tests [wavelet trees](../include/sdsl/wavelet_trees.hpp) on byte alphabets)
  * `wt-int-test` (tests [wavelet trees](../include/sdsl/wavelet_trees.hpp) on integer alphabets)

Test inputs are downloaded as needed before the first execution of the test.
See the [download.config](./download.config) files for details on the sources.

Executing `make test-sdsl` should take about 60 minutes on a recent machine.

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
