# Performance benchmark for sdsl data structures

This directory contains a set of benchmarks for [sdsl][sdsl]
data structures. Each benchmark is in its own subdirectory 
(e.g. `indexing_count`, `indexing_locate`) and
can be executed by changing in the specific directory and
typing `make timing`.

Directory `data` contains test inputs and can be populated
be executing the `get_corpus.sh` script. This will download
test inputs form the excellent [Pizza&Chili][pz] website.

Directory `tmp` is used to store temporary files (like
plain suffix arrays) which are used to generate compressed
structures.


## Author

Simon Gog (simon.gog@gmail.com)

[sdsl]: https://github.com/simongog/sdsl "sdsl"
[pz]: http://pizzachili.di.unipi.it "Pizza&Chili"
