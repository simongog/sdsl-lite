# Benchmarks for sdsl data structures

This directory contains a set of benchmarks for [sdsl][sdsl]
data structures. Each benchmark is in its own subdirectory and
so far we have:

* [indexing_count](./indexing_count): Evaluates the performance
  of count queries on different FM-Indexes/CSAs. Count query
  means _How many times occurs my pattern P in the text T?_
* [indexing_locate](./indexing_locate): Evaluates the performance
  of _locate queries_ on different FM-Indexes/CSAs. Locate query
  means _At which positions does pattern P occure in T?_

You can executed the benchmarks by going into benchmark's 
directory and typing `make timing`. But before doing this
you should populate the directory [data](./data) with test
inputs by executing the script 
[get_corpus.sh](./get_corpus.sh). This will download
test inputs form the excellent [Pizza&Chili][pz] website.

Directory `tmp` is used to store temporary files (like
plain suffix arrays) which are used to generate compressed
structures. 




## Literature

The benchmark code originates from the following article and can be used
to easily reproduce the results presented in the paper.


Simon Gog, Matthias Petri: _Optimized Succinct Data Structures for Massive Data_. 2013.
Accepted for publication in Software, Practice and Experience. 
[Preprint][http://people.eng.unimelb.edu.au/sgog/optimized.pdf]


## Author

Simon Gog (simon.gog@gmail.com)

[sdsl]: https://github.com/simongog/sdsl "sdsl"
[pz]: http://pizzachili.di.unipi.it "Pizza&Chili"
[PP]: http://people.eng.unimelb.edu.au/sgog/optimized.pdf "Preprint"
