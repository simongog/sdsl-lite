Directly Addressable Codes (DACs):
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

Directly Addressable Codes (DACs) consist in a variable-length encoding
scheme for integers that enables direct access to any element of the
encoded sequence and obtains compact spaces.

More information can be found in
http://lbd.udc.es/Repository/Publications/Drafts/DirAddVarLen.pdf
and
http://lbd.udc.es/Repository/Thesis/1303993824373_phdLadra.pdf

This version obtains the optimal space for the representation of
integers with no further restriction.

The software is under GNU/GPL v3 and comes without any warranty.



Source Code Included:
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

The folder src contains the source code for DACs.
More concretely, the use of DACs require the following files:
basics.c  basics.h  bitrankw32int.c  bitrankw32int.h  constants.h  dacs.c  dacs.h  

We also include a testing program, test.c, that generate a binary file, test,
to show how to use DACs' interface.


Compiling:
-=-=-=-=-

For compiling execute 'make'. With 'make clean' you can erase binary files.


Usage example:
-=-=-=-=-=-=-
For testing the source code use the binary file generated:
$ ./test <integer list file>  <outfile>


The input format for the integer list is a binary file containing a 4 bytes field
for each integer.

It saves DACs' representation of the integer list in <outfile>.
