# Tests for sdsl data structures

In directory contains test code for various data structures.
We used the gtest framework for testing. The test suite 
consists of test of basic to complex data structures.

You can compile the tests by calling `make` 
and execute them by calling `make test`
or `make vtest`, if you would like to
execute the tests with [valgrind][VG].



Test inputs  are downloaded as needed from a 
[Project Gutenberg][PG] mirror and Simon's
homepage before the first execution of the test.
Two tools have to be installed for the downloading
and extracting process:

 * [cURL][CURL] is required by the test input download script.
 * [gzip][GZIP] is required to extract compressed files.


Please report, if the test fail. Thanks. 

## Acknowledgements
  We thank 
  * Project Gutenberg for providing test inputs for the
    byte-alphabet data structures.
  * Shane Culpepper for providing the test inputs 
    `keeper.int` and `moby.int` for the integer-alphabet CSA and CSTs.

[VG]: http://valgrind.org/ "Valgrind"
[PG]: http://www.gutenberg.org/ "Project Gutenberg"
[CURL]: http://curl.haxx.se/ "cURL"
[GZIP]: http://www.gnu.org/software/gzip/ "Gzip Compressor"
