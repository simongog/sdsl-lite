INSTALL

1. git submodule update --init
2. mkdir build
3. cd build
4. cmake .. && make

CREATE COLLECTIONS

1. cd build
2. ./create-collection.x -c ../collection/NEWNAME -i NEWNAME.raw

EXECUTE BENCHMARK

1. cd build
2. ./gm_index-YOUR_IDX.x -c ../collections/your_collection
3. ./gm_search-YOUR_IDX.x -c ../collections/your_collection -p ../collections/your_collection/patterns/your_pattern.txt