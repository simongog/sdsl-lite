# Benchmarking operation `locate` on FM-indexes

## Methodology

Explored dimensions:
  
  * text type
  * instance size (just adjust the test_case.config file for this)
  * suffix array sampling density
  * index implementations

Pattern selection:

We use the methodology of [Ferragina et al.][FGNV08] (Section 5.3): 
,,Locate sufficient random patterns of length 5 to obtain a total of 2 to 
3 million occurrences''.

## Directory structure

  * [bin](./bin): Contains the executables of the project.
    * `build_idx_*` generates indexes
    * `query_idx_*` executes the locate experiments 
    * `info_*` outputs the space breakdown of an index.
    * `pattern_random` pattern generator.
  * [indexes](./indexes): Contains the generated indexes.
  * [pattern](./pattern): Contains the generated pattern.
  * [results](./results): Contains the results of the experiments.
  * [src](./src):  Contains the source code of the benchmark.
  * [visualize](./visualize): Contains a `R`-script which
			   generates a report in LaTeX format. 

	Files included in this archive from the [Pizza&Chili][pz] website:
      * [src/run_quries_sdsl.cpp](src/run_queries_sdsl.cpp) is a adapted version of the
	    Pizza&Chili file run_queries.c .

## Prerequisites
  * For the visualization you need the following software:
    - [R][RPJ] with package `tikzDevice`. You can install the
      package by calling 
      `install.packages("filehash", repos="http://cran.r-project.org")`
	  and 
	  `install.packages("tikzDevice", repos="http://R-Forge.R-project.org")`
	  in `R`.
    - Compressors [xz][XZ] and [gzip][GZIP] are used to get
	  compression baselines.
    - [pdflatex][LT] to generate the pdf reports.
  * The construction of the 200MB indexes requires about 1GB
    of RAM.
		
## Usage

 * `make timing` compiles the programs, downloads the 200MB
    [Pizza&Chili][pz] test cases, builds the indexes,
   runs the performance tests, and generated a report located at
   `visualize/locate.pdf`. The raw numbers of the timings
   can be found in the `results/all.txt`. 
   Indexes and temporary files are stored in the
   directory `indexes` and `tmp`. For the 5 x 200 MB of
   [Pizza&Chili][pz] data the project will produce about
   36 GB of additional data. On my machine (MacBookPro Retina
   2.6GHz Intel Core i7, 16GB 1600 Mhz DDR3, SSD) the
   benchmark, triggerd by `make timing`, took about 2 hours
   and 20 minutes (excluding the time to download the test instances).
   Have a look at the [generated report][RES].
 * All created indexes and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark
  The project contains several configuration files:
 
  * [index.config][IDXCONFIG]: Specify data structures' 
       ID, sdsl-class and LaTeX-name for the report.
  * [test_case.config][TCCONF]: Specify test cases's
       ID, path, LaTeX-name for the report, and download URL.
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
