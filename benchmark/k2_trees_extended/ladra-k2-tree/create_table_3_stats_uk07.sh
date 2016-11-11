        echo "Size in kBytes uncompressed"
        du -c table_3/uk07-05_2/*.il table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/uk07-05_2/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -c table_3/uk07-05_2/*.cil table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr table_3/uk07-05_2/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing uk07-05-bfs with 4 2 6 22"
        sleep 2
        tstime ./build_tree "/media/ssd/uk07-05-bfs.ladrabin" "table_3_bfs/uk07-05_2/uk07-05" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/uk07-05/*.il table_3_bfs/uk07-05/*.lv table_3_bfs/uk07-05_2/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/uk07-05/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/uk07-05/*.cil table_3_bfs/uk07-05_2/*.lv table_3_bfs/uk07-05_2/*.tr table_3_bfs/uk07-05_2/*.voc | tail -n 1 | awk '{print $1}'



        echo "Processing uk07-05-2-bfs"
        sleep 2
        tstime ./speed_direct table_3_bfs/uk07-05_2/uk07-05 /media/ssd/uk07-05.queries
                                                                                         
        echo "Processing uk07-05-2"
        sleep 2
        tstime ./speed_direct table_3/uk07-05_2/uk07-05 /media/ssd/uk07-05.queries

                echo "Processing uk07-05-2-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/uk07-05_2/uk07-05 /media/ssd/uk07-05.queries
                                                                                         
        echo "Processing uk07-05-2"
        sleep 2
        tstime ./speed_reveres table_3/uk07-05_2/uk07-05 /media/ssd/uk07-05.queries