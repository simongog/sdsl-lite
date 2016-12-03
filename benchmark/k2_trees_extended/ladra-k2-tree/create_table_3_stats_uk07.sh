
#!/bin/bash

print_garbage (){
  #fastest way to trick benchmark script
  echo "# comp_word_it = 0"
  echo "# frequency_encoding = 0"
  echo "# dac_compression = 0"
}

print_times (){
  #fastest way to trick benchmark script
  echo "# adj_time = 0"
  echo "# adj_check = 0"
  echo "# neighbors_time = 0"
  echo "# neighbors_check = 0"
  echo  "# reverse_neighbors_time = 0"
  echo "# reverse_neighbors_check = 0"
}

print_info (){
  echo "# K2_ID = LADRA_"$1
  echo "# TC_ID = "$1
  echo "# K2_TEX_NAME = LADRA_"$1
  echo "# TC_TEX_NAME = "$1
  echo "# TC_SIZE = 0"
}

  DATA_DIR="/media/ssd/"

  rm -rf table_3
  rm -rf table_3_bfs

	mkdir -p table_3/eu2005/
	mkdir -p table_3/indo/
	mkdir -p table_3/uk/
	mkdir -p table_3/arab/
	mkdir -p table_3/uk07-05/
	mkdir -p table_3/uk07-05_2/

  echo "Processing uk07-05 with 4 2 6 22"
  print_info "UK07-05"
  sleep 2
  ./build_tree ${DATA_DIR}"uk07-05.ladrabin" "table_3/uk07-05_2/uk07-05" 4 2 6 22 1
  du -cb table_3/uk07-05_2/*.il table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
  print_times
  ./compress_leaves "table_3/uk07-05_2/uk07-05" 20000000
  echo "Size in kBytes compressed"
  du -cb table_3/uk07-05_2/*.cil table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr table_3/uk07-05_2/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
  ./speed_direct table_3/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
  ./speed_reverse table_3/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
  print_garbage


  echo "Processing uk07-05-bfs with 4 2 6 22"
  sleep 2
  print_info "UK07-05-BFS"
  ./build_tree ${DATA_DIR}"uk07-05-bfs.ladrabin" "table_3_bfs/uk07-05_2/uk07-05" 4 2 6 22 1
  du -cb table_3_bfs/uk07-05/*.il table_3_bfs/uk07-05/*.lv table_3_bfs/uk07-05_2/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
  print_times
  ./compress_leaves "table_3_bfs/uk07-05/uk07-05" 20000000
  echo "Size in kBytes compressed"
  du -cb table_3_bfs/uk07-05/*.cil table_3_bfs/uk07-05_2/*.lv table_3_bfs/uk07-05_2/*.tr table_3_bfs/uk07-05_2/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
  ./speed_direct table_3_bfs/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
  ./speed_reverse table_3_bfs/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
  print_garbage
