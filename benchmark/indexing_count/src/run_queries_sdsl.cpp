/*
 * Run Queries
 */
#include <sdsl/suffix_arrays.hpp>
#include <string>

#include <stdlib.h>
#include "interface.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* only for getTime() */
#include <sys/time.h>
#include <sys/resource.h>

#define COUNT 		('C')
#define LOCATE 		('L')
#define EXTRACT 	('E')
#define DISPLAY 	('D')
#define VERBOSE 	('V')

using namespace sdsl;
using namespace std;

/* local headers */
void do_count(const CSA_TYPE&);
void do_locate(const CSA_TYPE&);
//void do_extract (void);
//void do_display(ulong length);
void pfile_info(ulong* length, ulong* numpatt);
//void output_char(uchar c, FILE * where);
double getTime(void);
void usage(char* progname);

static int Verbose = 0;
static ulong Index_size, Text_length;
static double Load_time;

/*
 * Temporary usage: run_queries <index file> <type> [length] [V]
 */
int main(int argc, char* argv[])
{
    char* filename;
    char querytype;

    if (argc < 2)	{
        usage(argv[0]);
        exit(1);
    }

    filename = argv[1];
    querytype = *argv[2];

    CSA_TYPE csa;
    fprintf(stderr, "# File = %s\n",(string(argv[1]) + "." + string(SUF)).c_str());
    fprintf(stderr, "# program = %s\n",string(SUF).c_str());
    Load_time = getTime();
    load_from_file(csa, (string(argv[1]) + "." + string(SUF)).c_str());
    Load_time = getTime() - Load_time;
    fprintf(stderr, "# Load_index_time_in_sec = %.2f\n", Load_time);
    std::cerr<<"# text_size = " << csa.size()-1 << std::endl;

    Index_size = size_in_bytes(csa);
    Text_length = csa.size()-1; // -1 since we added a sentinel character
    /*	Index_size /=1024; */
    fprintf(stderr, "# Index_size_in_bytes = %lu\n", Index_size);
    bool mapped = false;
#ifdef USE_HP
    mapped = mm::map_hp();
#endif
    fprintf(stderr, "# hugepages = %i\n", (int)mapped);
    bool use_sse = false;
#ifdef __SSE4_2__
    use_sse = true;
#endif
    fprintf(stderr, "# sse = %i\n", (int)use_sse);
    bool use_popcount_tl = false;
#ifdef POPCOUNT_TL
    use_popcount_tl = true;
#endif
    fprintf(stderr, "# popcount_tl = %i\n", (int)use_popcount_tl);

    switch (querytype) {
        case COUNT:
            if (argc > 3)
                if (*argv[3] == VERBOSE) {
                    Verbose = 1;
                    fprintf(stdout,"%c", COUNT);
                }
            do_count(csa);
            break;
        case LOCATE:
            if (argc > 3)
                if (*argv[3] == VERBOSE) {
                    Verbose = 1;
                    fprintf(stdout,"%c", LOCATE);
                }
            do_locate(csa);
            break;
            /*		case EXTRACT:
            			if (argc > 3)
            				if (*argv[3] == VERBOSE) {
            						Verbose = 1;
            						fprintf(stdout,"%c", EXTRACT);
            				}

            			do_extract();
            			break;
            		case DISPLAY:
            			if (argc < 4) {
            				usage(argv[0]);
            				exit (1);
            			}
            			if (argc > 4)
            				if (*argv[4] == VERBOSE){
            						Verbose = 1;
            						fprintf(stdout,"%c", DISPLAY);

            				}
            			do_display((ulong) atol(argv[3]));
            			break;
        */		default:
            fprintf(stderr, "Unknow option: main ru\n");
            exit(1);
    }
#ifdef USE_HP
    if (mapped) {
        mm::unmap_hp();
    }
#endif
    return 0;
}


void
do_count(const CSA_TYPE& csa)
{
    ulong numocc, length, tot_numocc = 0, numpatt, res_patt;
    double time, tot_time = 0;
    uchar* pattern;

    pfile_info(&length, &numpatt);
    res_patt = numpatt;

    pattern = (uchar*) malloc(sizeof(uchar) * (length));
    if (pattern == NULL) {
        fprintf(stderr, "Error: cannot allocate\n");
        exit(1);
    }

    while (res_patt) {

        if (fread(pattern, sizeof(*pattern), length, stdin) != length) {
            fprintf(stderr, "Error: cannot read patterns file\n");
            exit(1);
        }

        /* Count */
        time = getTime();
        numocc = sdsl::count(csa, pattern, pattern+length);

        if (Verbose) {
            fwrite(&length, sizeof(length), 1, stdout);
            fwrite(pattern, sizeof(*pattern), length, stdout);
            fwrite(&numocc, sizeof(numocc), 1, stdout);
        }
        tot_time += (getTime() - time);
        tot_numocc += numocc;
        res_patt--;
    }

    fprintf(stderr, "# Total_Num_occs_found = %lu\n", tot_numocc);
    fprintf(stderr, "# Count_time_in_milli_sec = %.4f\n", tot_time*1000);
    fprintf(stderr, "# Count_time/Pattern_chars = %.4f\n",
            (tot_time * 1000) / (length * numpatt));
    fprintf(stderr, "# Count_time/Num_patterns = %.4f\n\n",
            (tot_time * 1000) / numpatt);
    fprintf(stderr, "# (Load_time+Count_time)/Pattern_chars = %.4f\n",
            ((Load_time+tot_time) * 1000) / (length * numpatt));
    fprintf(stderr, "# (Load_time+Count_time)/Num_patterns = %.4f\n\n",
            ((Load_time+tot_time) * 1000) / numpatt);

    free(pattern);
}


void
do_locate(const CSA_TYPE& csa)
{
    ulong numocc, length; //, *occ,
    int_vector<32> occ;
    ulong tot_numocc = 0, numpatt;
    double time, tot_time = 0;
    uchar* pattern;

    pfile_info(&length, &numpatt);

    pattern = (uchar*) malloc(sizeof(uchar) * (length));
    if (pattern == NULL) {
        fprintf(stderr, "Error: cannot allocate\n");
        exit(1);
    }

    while (numpatt) {

        if (fread(pattern, sizeof(*pattern), length, stdin) != length) {
            fprintf(stderr, "Error: cannot read patterns file\n");
            perror("run_queries");
            exit(1);
        }
        // Locate
        time = getTime();
        numocc = sdsl::locate(csa, pattern, pattern+length, occ);
        tot_time += (getTime() - time);

        tot_numocc += numocc;
        numpatt--;

        if (Verbose) {
            fwrite(&length, sizeof(length), 1, stdout);
            fwrite(pattern, sizeof(*pattern), length, stdout);
            fwrite(&numocc, sizeof(numocc), 1, stdout);
        }

//		if (numocc) free (occ);
    }

    fprintf(stderr, "# Total_Num_occs_found = %lu\n", tot_numocc);
    fprintf(stderr, "# Locate_time_in_secs = %.2f\n", tot_time);
    fprintf(stderr, "# Locate_time/Num_occs = %.4f\n\n", (tot_time * 1000) / tot_numocc);
    fprintf(stderr, "# (Load_time+Locate_time)/Num_occs = %.4f\n\n", ((tot_time+Load_time) * 1000) / tot_numocc);

    free(pattern);
}


/*
void do_display(ulong numc) {

	int error = 0;
	ulong numocc, length, i, *snippet_len, tot_numcharext = 0, numpatt;
	double time, tot_time = 0;
	uchar *pattern, *snippet_text;

	pfile_info (&length, &numpatt);

	pattern = (uchar *) malloc (sizeof (uchar) * (length));
	if (pattern == NULL)
	{
		fprintf (stderr, "Error: cannot allocate\n");
		exit (1);
	}

	fprintf(stderr, "Snippet length %lu\n", numc);

	while (numpatt)
	{

		if (fread (pattern, sizeof (*pattern), length, stdin) != length)
		{
			fprintf (stderr, "Error: cannot read patterns file\n");
			perror ("run_queries");
			exit (1);
		}

		// Display
		time = getTime ();
		error =	display (Index, pattern, length, numc, &numocc,
				    	 &snippet_text, &snippet_len);
		IFERROR (error);
		tot_time += (getTime () - time);

		if (Verbose) {
			ulong j, len = length + 2*numc;
		    char blank = '\0';
			fwrite(&length, sizeof(length), 1, stdout);
			fwrite(pattern, sizeof(*pattern), length, stdout);
			fwrite(&numocc, sizeof(numocc), 1, stdout);
			fwrite(&len, sizeof(len), 1, stdout);

			for (i = 0; i < numocc; i++){
				fwrite(snippet_text+len*i,sizeof(uchar),snippet_len[i],stdout);
				for(j=snippet_len[i];j<len;j++)
				   fwrite(&blank,sizeof(uchar),1,stdout);
			}

		}
		numpatt--;

		for(i=0; i<numocc; i++) {
			tot_numcharext += snippet_len[i];
		}

		if (numocc) {
			free (snippet_len);
			free (snippet_text);
		}
	}

	fprintf (stderr, "#Total_num_chars_extracted = %lu\n", tot_numcharext);
	fprintf (stderr, "#Display_time_in_sec = %.2f secs\n", tot_time);
	fprintf (stderr, "#Time_display/Tot_num_chars = %.4f\n\n", (tot_time*1000) / tot_numcharext);
	fprintf (stderr, "#(Load_time+Time_display)/Tot_num_chars = %.4f\n\n", ((Load_time+tot_time)*1000) / tot_numcharext);

	free (pattern);
}
*/


/* Open patterns file and read header */
void
pfile_info(ulong* length, ulong* numpatt)
{
    int error;
    uchar c;
    uchar origfilename[257];

    error = fscanf(stdin, "# number=%lu length=%lu file=%s forbidden=", numpatt,
                   length, origfilename);
    if (error != 3) {
        fprintf(stderr, "Error: Patterns file header not correct\n");
        perror("run_queries");
        exit(1);
    }

    fprintf(stderr, "# pat_cnt = %lu\n", *numpatt);
    fprintf(stderr, "# pat_length = %lu\n", *length);
    fprintf(stderr, "# forbidden_chars = ");

    while ((c = fgetc(stdin)) != 0) {
        if (c == '\n') break;
        fprintf(stderr, "%d",c);
    }

    fprintf(stderr, "\n");

}

/*
void
do_extract()
{
	int error = 0;
	uchar *text, orig_file[257];
	ulong num_pos, from, to, numchars, readen, tot_ext = 0;
	double time, tot_time = 0;

	error = fscanf(stdin, "# number=%lu length=%lu file=%s\n", &num_pos, &numchars, orig_file);
	if (error != 3)
	{
		fprintf (stderr, "Error: Intervals file header is not correct\n");
		perror ("run_queries");
		exit (1);
	}
	fprintf(stderr, "# number=%lu length=%lu file=%s\n", num_pos, numchars, orig_file);

	while(num_pos) {

		if(fscanf(stdin,"%lu,%lu\n", &from, &to) != 2) {
			fprintf(stderr, "Cannot read correctly intervals file\n");
			exit(1);
		}

		time = getTime();
		error = extract(Index, from, to, &text, &readen);
		IFERROR(error);
		tot_time += (getTime() - time);

		tot_ext += readen;

		if (Verbose) {
			fwrite(&from,sizeof(ulong),1,stdout);
			fwrite(&readen,sizeof(ulong),1,stdout);
			fwrite(text,sizeof(uchar),readen, stdout);
		}

		num_pos--;
		free(text);
	}

	fprintf (stderr, "# Total_num_chars_extracted = %lu\n", tot_ext);
	fprintf (stderr, "# Extract_time_in_sec = %.2f\n", tot_time);
	fprintf (stderr, "# Extract_time/Num_chars_extracted = %.4f\n\n",
		 (tot_time * 1000) / tot_ext);
	fprintf (stderr, "(Load_time+Extract_time)/Num_chars_extracted = %.4f\n\n",
		 ((Load_time+tot_time) * 1000) / tot_ext);
}
*/

double
getTime(void)
{

    double usertime, systime;
    struct rusage usage;

    getrusage(RUSAGE_SELF, &usage);

    usertime = (double) usage.ru_utime.tv_sec +
               (double) usage.ru_utime.tv_usec / 1000000.0;

    systime = (double) usage.ru_stime.tv_sec +
              (double) usage.ru_stime.tv_usec / 1000000.0;

    return (usertime + systime);

}

void usage(char* progname)
{
    fprintf(stderr, "\nThe program loads <index> and then executes over it the\n");
    fprintf(stderr, "queries it receives from the standard input. The standard\n");
    fprintf(stderr, "input comes in the format of the files written by \n");
    fprintf(stderr, "genpatterns or genintervals.\n");
    fprintf(stderr, "%s reports on the standard error time statistics\n", progname);
    fprintf(stderr, "regarding to running the queries.\n\n");
    fprintf(stderr, "Usage:  %s <index> <type> [length] [V]\n", progname);
    fprintf(stderr, "\n\t<type>   denotes the type of queries:\n");
    fprintf(stderr, "\t         %c counting queries;\n", COUNT);
    fprintf(stderr, "\t         %c locating queries;\n", LOCATE);
    fprintf(stderr, "\t         %c displaying queries;\n", DISPLAY);
    fprintf(stderr, "\t         %c extracting queries.\n\n", EXTRACT);
    fprintf(stderr, "\n\t[length] must be provided in case of displaying queries (D)\n");
    fprintf(stderr, "\t         and denotes the number of characters to display\n");
    fprintf(stderr, "\t         before and after each pattern occurrence.\n");
    fprintf(stderr, "\n\t[V]      with this options it reports on the standard output\n");
    fprintf(stderr, "\t         the results of the queries. The results file should be\n");
    fprintf(stderr, "\t         compared with trusted one by compare program.\n\n");
}
