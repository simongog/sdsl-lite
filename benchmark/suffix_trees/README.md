# Benchmarking suffix trees

## Methodology

Explored dimensions:
  
  * suffix tree implementations
  * test cases
  * methods (`LCA`, `Letter`, `SLink`, `Child`, `SDepth`, `Parent`)

## Directory structure

  * [bin](./bin): Contains the executables of the project.
  * [indexes](./indexes): Contains the serialized suffix trees.
  * [results](./results): Contains the results of the experiments.
  * [src](./src): Contains the source code of the benchmark.
  * [stats](./stats): Contains statistics about the construction and
					  the structure of the indexes.
  * [visualize](./visualize): Contains a `R`-script which generates
							  a report in LaTeX format.

## Prerequisites

  * For the visualization you need the following software:
    - [R][RPJ] to generate the tex files.
    - [pdflatex][LT] to generate the pdf reports.

## Usage

 * `make timing` compiles the programs, downloads the test instances,
   builds the suffix trees, runs the performance tests and generates a
   report located at `visualize/suffix_trees.pdf`. The raw numbers of
   the timings can be found in the `results/all.txt` file.
 * All created binaries and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark

The project contains several configuration files:
 
  * [index.config][STCONFIG]: Specify different suffix tree implementations.
  * [test_case.config][TCCONF]: Specify test instances by ID, path, LaTeX-name 
								for the report, and download URL.
  * [compile_options.config][CCONF]: Specify compile options by option string.

Note that the benchmark will execute every combination of suffix trees and test cases.

[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[STCONFIG]: ./index.config "index.config"
[TCCONF]: ./test_case.config "test_case.config"
[CCONF]: ./compile_options.config "compile_options.config"

