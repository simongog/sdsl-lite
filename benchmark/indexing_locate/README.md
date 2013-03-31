# Benchmarking operation `locate` on FM-indexes



## Methodology

Explored Dimensions:
  
  * text type
  * text size (just adjust the test_case.config file for this)
  * suffix array sampling density
  * index implementations

Pattern selection:

We use the methodology of [Ferragina et al.][FGNV08] (Section 5.3): 
,,Locate sufficient random patterns of length 5 to obtain a total of 2 to 
3 million occurrences''.

## Directory structure

  * [bin](./bin): Contains the executables of the project.
    * `build_idx_*` generates indexes
    * `query_idx_*` executes the count experiments 
    * `info_*` outputs the space breakdown of an index.
    * `pattern_random` pattern generator.
  * [src](./src):  Contains the source code of the benchmark.
  * [visualize](./visualize): `R` and `pdflatex` scripts to generate
               reports of the collected data.
  * [results](./results): Contains the results of the experiments.

	Files included in this archive form the Pizza&Chili website:
      * [src/run_quries_sdsl.cpp](src/run_queries_sdsl.cpp) is a adapted version of the
	    Pizza&Chili file run_queries.c .

## Prerequisites
  * For the visualization you need the following software:
    - [R][RPJ] with package `tikzdevice`. You can install the
      package by calling `install.packages("tikzdevice")` in R.
    - Compressors [xz][XZ] and [gzip][GZIP] are used to get
	  compression baselines.
    - [pdflatex][LT] to generate the pdf reports.
  * The construction of the 200MB indexes needs about 1GB
    of RAM.
		
## Usage

 * Populate the data directory with some test files. E.g. you
   can download the 200MB [Pizza&Chili][pz] test cases by 
   executing [../data/get_corpus.sh](../data/get_corpus.sh). 
 * `make timing`  compiles the programs, builds the indexes and
   runs the performance tests and generated a report located as
   `visualization/count.pdf`. The raw numbers of the timing
   can be found in the `results/all.txt`. 
   Indexes and temporary files are stored in the
   directory `data` and `tmp`. For the 5 x 200 MB of
   [Pizza&Chili][pz] data the project will produce about
   36 GB of additional data. On my machine (MacBookPro Retina
   2.6GHz Intel Core i7, 16GB 1600 Mhz DDR3, SSD) the
   benchmark, invoced by `make timing`, took about 2 hours.
   Have a look at the [generated report][RES].
 * All created indexes and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark
  The project contains several configuration files:
 
  * [index.config][IDXCONFIG]: Specify data structures' 
       ID, sdsl-class and LaTeX-name for the report.
  * [test_case.config][TCCONF]: Specify test cases's
       ID, path, and LaTeX-name for the report.
  * [sample.config][SCONF]: Specify samplings' ID,
       rate for SA, and rate for ISA. 

  Note that the benchmark will execute every combination of your
  choices. 

  Finally, the visualization can also be configured:

  * [visualize/index-filter.config][VCONF]: Specify which 
   indexes should be listed in the report and which style should be used.

[sdsl]: https://github.com/simongog/sdsl "sdsl"
[pz]: http://pizzachili.di.unipi.it "Pizza&Chili"
[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[RES]: https://github.com/simongog/simongog.github.com/raw/master/assets/images/locate.pdf "locate.pdf"
[FGNV08]: http://dl.acm.org/citation.cfm?doid=1412228.1455268 "FGNV08"
[IDXCONFIG]: ./index.config "index.config"
[TCCONF]: ./test_case.config "test_case.config"
[SCONF]: ./sample.config "sample.config"
[VCONF]: ./visualize/index-filter.config "index-filter.config"
[XZ]: http://tukaani.org/xz/ "XZ Compressor"
[GZIP]: http://www.gnu.org/software/gzip/ "Gzip Compressor"
