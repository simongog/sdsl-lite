
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

        DATA_DIR="/home/d056848/webgraphs/"

    rm -rf table_3
    rm -rf table_3_bfs

	mkdir -p table_3/eu2005/
	mkdir -p table_3/indo/
	mkdir -p table_3/uk/
	mkdir -p table_3/arab/
	mkdir -p table_3/uk07-05/
	mkdir -p table_3/uk07-05_2/

	mkdir -p table_3_bfs/eu2005/
	mkdir -p table_3_bfs/indo/
	mkdir -p table_3_bfs/uk/
	mkdir -p table_3_bfs/arab/
	mkdir -p table_3_bfs/uk07-05/
	mkdir -p table_3_bfs/uk07-05_2/
        echo "Processing eu2005"
        sleep 2
        print_info "EU2005"
        ./build_tree ${DATA_DIR}"eu2005.ladrabin" "table_3/eu2005/eu2005" 4 2 5 18 1
        du -cb table_3/eu2005/*.il table_3/eu2005/*.lv table_3/eu2005/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3/eu2005/eu2005" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3/eu2005/*.cil table_3/eu2005/*.lv table_3/eu2005/*.tr table_3/eu2005/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3/eu2005/eu2005 ${DATA_DIR}eu2005.queries.single 1
        ./speed_direct table_3/eu2005/eu2005 ${DATA_DIR}eu2005.queries 1
        ./speed_reverse table_3/eu2005/eu2005 ${DATA_DIR}eu2005.queries 1
        print_garbage


        echo "Processing indo"
        print_info "INDO"
        sleep 2
        ./build_tree ${DATA_DIR}"indo.ladrabin" "table_3/indo/indo" 4 2 5 20 1
        du -cb table_3/indo/*.il table_3/indo/*.lv table_3/indo/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3/indo/indo" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3/indo/*.cil table_3/indo/*.lv table_3/indo/*.tr table_3/indo/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3/indo/indo ${DATA_DIR}indo.queries.single 1
        ./speed_direct table_3/indo/indo ${DATA_DIR}indo.queries 1
        ./speed_reverse table_3/indo/indo ${DATA_DIR}indo.queries 1
        print_garbage

        echo "Processing uk"
        print_info "UK2002"
        sleep 2
        ./build_tree ${DATA_DIR}"uk.ladrabin" "table_3/uk/uk" 4 2 6 22 1
        du -cb table_3/uk/*.il table_3/uk/*.lv table_3/uk/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3/uk/uk" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3/uk/*.cil table_3/uk/*.lv table_3/uk/*.tr table_3/uk/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3/uk/uk ${DATA_DIR}uk.queries.single 1
        ./speed_direct table_3/uk/uk ${DATA_DIR}uk.queries 1
        ./speed_reverse table_3/uk/uk ${DATA_DIR}uk.queries 1
        print_garbage
                
        echo "Processing arab"
        print_info "ARAB"
        sleep 2
        ./build_tree ${DATA_DIR}"arab.ladrabin" "table_3/arab/arab" 4 2 6 22 1
        du -cb table_3/arab/*.il table_3/arab/*.lv table_3/arab/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3/arab/arab" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3/arab/*.cil table_3/arab/*.lv table_3/arab/*.tr table_3/arab/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3/arab/arab ${DATA_DIR}arab.queries.single 1
        ./speed_direct table_3/arab/arab ${DATA_DIR}arab.queries 1
        ./speed_reverse table_3/arab/arab ${DATA_DIR}arab.queries 1
        print_garbage

        echo "Processing uk07-05 with 4 2 6 22"
        print_info "UK07-05"
        sleep 2
        ./build_tree ${DATA_DIR}"uk07-05.ladrabin" "table_3/uk07-05_2/uk07-05" 4 2 6 22 1
        du -cb table_3/uk07-05_2/*.il table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3/uk07-05_2/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -cb table_3/uk07-05_2/*.cil table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr table_3/uk07-05_2/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries.single 1
        ./speed_direct table_3/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
        ./speed_reverse table_3/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
        print_garbage


        echo "Processing eu2005-bfs"
        sleep 2
        print_info "EU2005-BFS"
        ./build_tree ${DATA_DIR}"eu2005-bfs.ladrabin" "table_3_bfs/eu2005/eu2005" 4 2 5 18 1
        du -cb table_3_bfs/eu2005/*.il table_3_bfs/eu2005/*.lv table_3_bfs/eu2005/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3_bfs/eu2005/eu2005" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3_bfs/eu2005/*.cil table_3_bfs/eu2005/*.lv table_3_bfs/eu2005/*.tr table_3_bfs/eu2005/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3_bfs/eu2005/eu2005 ${DATA_DIR}eu2005.queries.single 1
        ./speed_direct table_3_bfs/eu2005/eu2005 ${DATA_DIR}eu2005.queries 1
        ./speed_reverse table_3_bfs/eu2005/eu2005 ${DATA_DIR}eu2005.queries 1
        print_garbage

        echo "Processing indo-bfs"
        sleep 2
        print_info "INDO-BFS"
        ./build_tree ${DATA_DIR}"indo-bfs.ladrabin" "table_3_bfs/indo/indo" 4 2 5 20 1
        du -cb table_3_bfs/indo/*.il table_3_bfs/indo/*.lv table_3_bfs/indo/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3_bfs/indo/indo" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3_bfs/indo/*.cil table_3_bfs/indo/*.lv table_3_bfs/indo/*.tr table_3_bfs/indo/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3_bfs/indo/indo ${DATA_DIR}indo.queries.single 1
        ./speed_direct table_3_bfs/indo/indo ${DATA_DIR}indo.queries 1
        ./speed_reverse table_3_bfs/indo/indo ${DATA_DIR}indo.queries 1
        print_garbage

        echo "Processing uk-bfs"
        sleep 2
        print_info "UK2002-BFS"
        ./build_tree ${DATA_DIR}"uk-bfs.ladrabin" "table_3_bfs/uk/uk" 4 2 6 22 1
        du -cb table_3_bfs/uk/*.il table_3_bfs/uk/*.lv table_3_bfs/uk/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3_bfs/uk/uk" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3_bfs/uk/*.cil table_3_bfs/uk/*.lv table_3_bfs/uk/*.tr table_3_bfs/uk/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3_bfs/uk/uk ${DATA_DIR}uk.queries.single 1
        ./speed_direct table_3_bfs/uk/uk ${DATA_DIR}uk.queries 1
        ./speed_reverse table_3_bfs/uk/uk ${DATA_DIR}uk.queries 1
        print_garbage
                
        echo "Processing arab-bfs"
        sleep 2
        print_info "ARAB-BFS"
        ./build_tree ${DATA_DIR}"arab-bfs.ladrabin" "table_3_bfs/arab/arab" 4 2 6 22 1
        du -cb table_3_bfs/arab/*.il table_3_bfs/arab/*.lv table_3_bfs/arab/*.tr | tail -n 1 | awk '{print "# uncompressed_size = "$1}'
        print_times
        ./compress_leaves "table_3_bfs/arab/arab" 4000000
        echo "Size in kBytes compressed"
        du -cb table_3_bfs/arab/*.cil table_3_bfs/arab/*.lv table_3_bfs/arab/*.tr table_3_bfs/arab/*.voc | tail -n 1 | awk '{print "# compressed_size = "$1}'
        ./speed_link table_3_bfs/arab/arab ${DATA_DIR}arab.queries.single 1
        ./speed_direct table_3_bfs/arab/arab ${DATA_DIR}arab.queries 1
        ./speed_reverse table_3_bfs/arab/arab ${DATA_DIR}arab.queries 1
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
        ./speed_link table_3_bfs/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries.single 1
        ./speed_direct table_3_bfs/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
        ./speed_reverse table_3_bfs/uk07-05_2/uk07-05 ${DATA_DIR}uk07-05.queries 1
        print_garbage
