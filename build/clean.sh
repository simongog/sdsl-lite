#!/bin/bash
# activate globbing
shopt -s extglob
# removes all but the listed files in the directory
# and subdirectories
rm -r !(build.sh|clean.sh|.gitignore)  

