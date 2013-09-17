#!/bin/bash
# activate globbing
shopt -s extglob
# removes all but the listed files in the directory
# and subdirectories
rm -rf !(run.sh|clean.sh|.gitignore|README.md)  

