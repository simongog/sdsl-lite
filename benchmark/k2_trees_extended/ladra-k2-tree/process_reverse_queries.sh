	#!/bin/bash
        echo "Processing eu2005"
        sleep 2
        tstime ./speed_reverse table_3/eu2005/eu2005 /media/ssd/eu2005.queries

        echo "Processing indo"
        sleep 2
        tstime ./speed_reverse table_3/indo/indo /media/ssd/indo.queries

        echo "Processing uk"
        sleep 2
        tstime ./speed_reverse table_3/uk/uk /media/ssd/uk.queries

        echo "Processing arab"
        sleep 2
        tstime ./speed_reverse table_3/arab/arab /media/ssd/arab.queries

        echo "Processing uk07-05"
        sleep 2
        tstime ./speed_reverse table_3/uk07-05/uk07-05 /media/ssd/uk07-05.queries

        echo "Processing uk07-05-2"
        sleep 2
        tstime ./speed_reverse table_3/uk07-05_2/uk07-05 /media/ssd/uk07-05.queries

        echo "Processing eu2005-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/eu2005/eu2005 /media/ssd/eu2005.queries

        echo "Processing indo-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/indo/indo /media/ssd/indo.queries

        echo "Processing uk-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/uk/uk /media/ssd/uk.queries

        echo "Processing arab-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/arab/arab /media/ssd/arab.queries

        echo "Processing uk07-05-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/uk07-05/uk07-05 /media/ssd/uk07-05.queries

        echo "Processing uk07-05-2-bfs"
        sleep 2
        tstime ./speed_reverse table_3_bfs/uk07-05_2/uk07-05 /media/ssd/uk07-05.queries