#!/usr/bin/python

import fnmatch
import sys
import os
import re

includes = ['*.cpp','*.hpp'] # sources only
includes = r'|'.join([fnmatch.translate(x) for x in includes])

collection_path = "collection.txt"

def main():

    if len(sys.argv) == 1:
        print "Usage ./", sys.argv[0], "directory"
        sys.exit(0)
    cur_dir = sys.argv[1]
    doc_paths = []
    for root, dirs, files in os.walk(cur_dir):
        files = [f for f in files if re.match(includes, f)]
        for f in files:
            doc_path = root+"/"+f
            doc_paths.append(doc_path)

    print "Found ", len(doc_paths), "source files in", cur_dir
    collection_f = open(collection_path,'w')
    for doc_path in doc_paths:
        f = open(doc_path, 'r')
        content = f.read()
        collection_f.write(content)
        collection_f.write('\1')

if __name__ == '__main__':
    main()
