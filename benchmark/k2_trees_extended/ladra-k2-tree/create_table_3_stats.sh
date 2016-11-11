#run in su bash!!!!
#!/bin/bash
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
        tstime ./build_tree "/media/ssd/eu2005.ladrabin" "table_3/eu2005/eu2005" 4 2 5 18 1
        echo "Size in kBytes uncompressed"
        du -c table_3/eu2005/*.il table_3/eu2005/*.lv table_3/eu2005/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/eu2005/eu2005" 4000000
        echo "Size in kBytes compressed"
        du -c table_3/eu2005/*.cil table_3/eu2005/*.lv table_3/eu2005/*.tr table_3/eu2005/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing indo"
        sleep 2
        tstime ./build_tree "/media/ssd/indo.ladrabin" "table_3/indo/indo" 4 2 5 20 1
        echo "Size in kBytes uncompressed"
        du -c table_3/indo/*.il table_3/indo/*.lv table_3/indo/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/indo/indo" 4000000
        echo "Size in kBytes compressed"
        du -c table_3/indo/*.cil table_3/indo/*.lv table_3/indo/*.tr table_3/indo/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing uk"
        sleep 2
        tstime ./build_tree "/media/ssd/uk.ladrabin" "table_3/uk/uk" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3/uk/*.il table_3/uk/*.lv table_3/uk/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/uk/uk" 4000000
        echo "Size in kBytes compressed"
        du -c table_3/uk/*.cil table_3/uk/*.lv table_3/uk/*.tr table_3/uk/*.voc | tail -n 1 | awk '{print $1}'
		
	echo "Processing arab"
        sleep 2
        tstime ./build_tree "/media/ssd/arab.ladrabin" "table_3/arab/arab" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3/arab/*.il table_3/arab/*.lv table_3/arab/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/arab/arab" 4000000
        echo "Size in kBytes compressed"
        du -c table_3/arab/*.cil table_3/arab/*.lv table_3/arab/*.tr table_3/arab/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing uk07-05"
        sleep 2
        tstime ./build_tree "/media/ssd/uk07-05.ladrabin" "table_3/uk07-05/uk07-05" 4 2 8 24 1
        echo "Size in kBytes uncompressed"
        du -c table_3/uk07-05/*.il table_3/uk07-05/*.lv table_3/uk07-05/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/uk07-05/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -c table_3/uk07-05/*.cil table_3/uk07-05/*.lv table_3/uk07-05/*.tr table_3/uk07-05/*.voc | tail -n 1 | awk '{print $1}'

	echo "Processing uk07-05 with 4 2 6 22"
        sleep 2
        tstime ./build_tree "/media/ssd/uk07-05.ladrabin" "table_3/uk07-05_2/uk07-05" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3/uk07-05_2/*.il table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3/uk07-05_2/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -c table_3/uk07-05_2/*.cil table_3/uk07-05_2/*.lv table_3/uk07-05_2/*.tr table_3/uk07-05_2/*.voc | tail -n 1 | awk '{print $1}'


                echo "Processing eu2005-bfs"
        sleep 2
        tstime ./build_tree "/media/ssd/eu2005-bfs.ladrabin" "table_3_bfs/eu2005/eu2005" 4 2 5 18 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/eu2005/*.il table_3_bfs/eu2005/*.lv table_3_bfs/eu2005/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/eu2005/eu2005" 4000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/eu2005/*.cil table_3_bfs/eu2005/*.lv table_3_bfs/eu2005/*.tr table_3_bfs/eu2005/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing indo-bfs"
        sleep 2
        tstime ./build_tree "/media/ssd/indo-bfs.ladrabin" "table_3_bfs/indo/indo" 4 2 5 20 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/indo/*.il table_3_bfs/indo/*.lv table_3_bfs/indo/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/indo/indo" 4000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/indo/*.cil table_3_bfs/indo/*.lv table_3_bfs/indo/*.tr table_3_bfs/indo/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing uk-bfs"
        sleep 2
        tstime ./build_tree "/media/ssd/uk-bfs.ladrabin" "table_3_bfs/uk/uk" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/uk/*.il table_3_bfs/uk/*.lv table_3_bfs/uk/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/uk/uk" 4000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/uk/*.cil table_3_bfs/uk/*.lv table_3_bfs/uk/*.tr table_3_bfs/uk/*.voc | tail -n 1 | awk '{print $1}'
		
	echo "Processing arab-bfs"
        sleep 2
        tstime ./build_tree "/media/ssd/arab-bfs.ladrabin" "table_3_bfs/arab/arab" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/arab/*.il table_3_bfs/arab/*.lv table_3_bfs/arab/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/arab/arab" 4000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/arab/*.cil table_3_bfs/arab/*.lv table_3_bfs/arab/*.tr table_3_bfs/arab/*.voc | tail -n 1 | awk '{print $1}'

        echo "Processing uk07-05-bfs"
        sleep 2
        tstime ./build_tree "/media/ssd/uk07-05-bfs.ladrabin" "table_3_bfs/uk07-05/uk07-05" 4 2 8 24 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/uk07-05/*.il table_3_bfs/uk07-05/*.lv table_3_bfs/uk07-05/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/uk07-05/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/uk07-05/*.cil table_3_bfs/uk07-05/*.lv table_3_bfs/uk07-05/*.tr table_3_bfs/uk07-05/*.voc | tail -n 1 | awk '{print $1}'

	echo "Processing uk07-05-bfs with 4 2 6 22"
        sleep 2
        tstime ./build_tree "/media/ssd/uk07-05-bfs.ladrabin" "table_3_bfs/uk07-05_2/uk07-05" 4 2 6 22 1
        echo "Size in kBytes uncompressed"
        du -c table_3_bfs/uk07-05/*.il table_3_bfs/uk07-05/*.lv table_3_bfs/uk07-05_2/*.tr | tail -n 1 | awk '{print $1}'
        tstime ./compress_leaves "table_3_bfs/uk07-05/uk07-05" 20000000
        echo "Size in kBytes compressed"
        du -c table_3_bfs/uk07-05/*.cil table_3_bfs/uk07-05_2/*.lv table_3_bfs/uk07-05_2/*.tr table_3_bfs/uk07-05_2/*.voc | tail -n 1 | awk '{print $1}'
