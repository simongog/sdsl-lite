# Visualize the output of the resource logger

This directory contains a shell script [mm-log.sh](./mm-log.sh)
which can be used to visualize the output of the resource 
logger. See [mm-log.cpp](../../tutorial/mm-log.cpp) for
an example program which uses the logger.
An example output produced with this program can be found
in [english-200MB.mem-log](../english-200MB.mem-log) 
(the log for the construction of a CST for the Pizza&Chili 
 file `english.200MB`).

A call of `mm-log.sh english-200MB.mem-log` will call the
R-script `mm-log.R`, which will produce a resource usage
plot. The plot is stored in `mm-log.pdf`.

The coloring of the different phases of the CST construction
can be configured by editing the [mm-log.config](./mm-log.config) file.


