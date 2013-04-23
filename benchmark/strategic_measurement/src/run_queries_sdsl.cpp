/*
 * Run Queries
 */
#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <algorithm>

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

/* macro to detect and to notify errors */
#define IFERROR(error) {{if (error) { fprintf(stderr, "%s\n", error_index(error)); exit(1); }}}

using namespace sdsl;
using namespace std;

/* local headers */
void do_count(const CSA_TYPE&);
void do_locate(const CSA_TYPE&);
void do_extract(const CSA_TYPE&);
//void do_display(ulong length);
void pfile_info(ulong* length, ulong* numpatt);
//void output_char(uchar c, FILE * where);
double getTime(void);
void usage(char* progname);

static void* Index;	 /* opaque data type */
static int Verbose = 0;
static ulong Index_size, Text_length;
static double Load_time;
static ulong MAX_PAT=100000000;
static ulong SKIP=0;

/*
 * Temporary usage: run_queries <index file> <type> [length] [V]
 */
int main(int argc, char* argv[])
{
    int error = 0;
    char* filename;
    char querytype;

    if (argc < 2)	{
        usage(argv[0]);
        exit(1);
    }

    filename = argv[1];
    querytype = *argv[2];

    CSA_TYPE csa;
    fprintf(stderr, "Load from file %s\n",argv[1]);
    Load_time = getTime();
    load_from_file(csa, string(argv[1]));
    IFERROR(error);
    Load_time = getTime() - Load_time;
    fprintf(stderr, "# Load_index_time_in_sec = %.2f\n", Load_time);
    fprintf(stderr, "# text_size = %lld\n", csa.size()-1);

    Index_size = size_in_bytes(csa);
    IFERROR(error);
    Text_length = csa.size()-1; // -1 since we added a sentinel character
//	error = get_length(Index, &Text_length);
    IFERROR(error);
    /*	Index_size /=1024; */
    fprintf(stderr, "# Index_size_in_bytes = %lu\n", Index_size);
#ifdef USE_HP
    bool mapped = mm::map_hp();
    fprintf(stderr, "# hugepages = %i\n", (int)mapped);
#endif

    switch (querytype) {
        case COUNT:
            if (argc > 3) {
                MAX_PAT = atoi(argv[3]);
            }
            if (argc > 4) {
                SKIP = atoi(argv[4]);
            }
            do_count(csa);
            break;
        case LOCATE:
            if (argc > 3) {
                MAX_PAT = atoi(argv[3]);
            }
            if (argc > 4) {
                SKIP = atoi(argv[4]);
            }
            do_locate(csa);
            break;
        case EXTRACT:
            do_extract(csa);
            break;
        default:
            fprintf(stderr, "Unknow option: main ru\n");
            exit(1);
    }
#ifdef USE_HP
    if (mapped) {
        mm::unmap_hp();
    }
#endif

//	error = free_index(Index);
    IFERROR(error);

    return 0;
}


void
do_count(const CSA_TYPE& csa)
{
    int error = 0;
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
            perror("run_queries");
            exit(1);
        }

        /* Count */
        time = getTime();
        numocc = count(csa, pattern, pattern+length);
        IFERROR(error);

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
    int error = 0;
    ulong numocc, length; //, *occ,
    int_vector<32> occ;
    ulong tot_numocc = 0, numpatt = 0, processed_pat = 0;
    double time, tot_time = 0;
    uchar* pattern;

    pfile_info(&length, &numpatt);

    pattern = (uchar*) malloc(sizeof(uchar) * (length));
    if (pattern == NULL) {
        fprintf(stderr, "Error: cannot allocate\n");
        exit(1);
    }
    /*SG: added timeout of 60 seconds */
//    while (numpatt and tot_time < 60.0) {
    while (numpatt) {

        if (fread(pattern, sizeof(*pattern), length, stdin) != length) {
            fprintf(stderr, "Error: cannot read patterns file\n");
            perror("run_queries");
            exit(1);
        }
        // Locate
        time = getTime();
        numocc = locate(csa, pattern, pattern+length, occ);
        IFERROR(error);
        tot_time += (getTime() - time);
        ++processed_pat;

        tot_numocc += numocc;
        numpatt--;

        if (Verbose) {
            fwrite(&length, sizeof(length), 1, stdout);
            fwrite(pattern, sizeof(*pattern), length, stdout);
            fwrite(&numocc, sizeof(numocc), 1, stdout);
        }
    }

    fprintf(stderr, "# processed_pattern = %lu\n", processed_pat);
    fprintf(stderr, "# Total_Num_occs_found = %lu\n", tot_numocc);
    fprintf(stderr, "# Locate_time_in_secs = %.5f\n", tot_time);
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

    error = fscanf(stdin, "# number=%lu length=%lu file=%s ", numpatt,
                   length, origfilename);
    ulong numpatt2 = std::min(*numpatt, MAX_PAT);
    if (error != 3) {
        fprintf(stderr, "Error: Patterns file header not correct\n");
        perror("run_queries");
        exit(1);
    }

    fprintf(stderr, "# file_pat_cnt = %lu\n", *numpatt);
    fprintf(stderr, "# pat_cnt = %lu\n", numpatt2);
    fprintf(stderr, "# pat_length = %lu\n", *length);

    while ((c = fgetc(stdin)) != 0) {
        if (c == '\n') break;
    }
    fprintf(stderr, "\n");
    if (SKIP > 0 and SKIP + numpatt2 <= *numpatt) {
        uchar* pattern = new uchar[(*length)*SKIP];
        if (fread(pattern, sizeof(uchar), (*length)*SKIP, stdin) != SKIP*(*length)) {
            perror("could not skip patterns");
            delete [] pattern;
            exit(1);
        }
        delete [] pattern;
        fprintf(stderr, "# skip = %lu\n", SKIP);
    } else {
        fprintf(stderr, "# skip = %lu\n", 0);
    }
    *numpatt = numpatt2;
}

void
do_extract(const CSA_TYPE& csa)
{
    int error = 0;
    uchar* text, orig_file[257];
    ulong num_pos, from, to, numchars, tot_ext = 0;
    CSA_TYPE::size_type readlen = 0;
    double time, tot_time = 0;

    error = fscanf(stdin, "# number=%lu length=%lu file=%s\n", &num_pos, &numchars, orig_file);
    if (error != 3) {
        fprintf(stderr, "Error: Intervals file header is not correct\n");
        perror("run_queries");
        exit(1);
    }
    fprintf(stderr, "# number=%lu length=%lu file=%s\n", num_pos, numchars, orig_file);

    while (num_pos) {

        if (fscanf(stdin,"%lu,%lu\n", &from, &to) != 2) {
            fprintf(stderr, "Cannot read correctly intervals file\n");
            exit(1);
        }

        time = getTime();
        text = (uchar*)malloc(to-from+2);
        readlen = sdsl::extract(csa, from, to, text);
        tot_time += (getTime() - time);

        tot_ext += readlen;

        if (Verbose) {
            fwrite(&from,sizeof(ulong),1,stdout);
            fwrite(&readlen,sizeof(ulong),1,stdout);
            fwrite(text,sizeof(uchar),readlen, stdout);
        }

        readlen = 0;
        num_pos--;
        free(text);
    }

    fprintf(stderr, "# Total_num_chars_extracted = %lu\n", tot_ext);
    fprintf(stderr, "# Extract_time_in_sec = %.2f\n", tot_time);
    fprintf(stderr, "# Extract_time/Num_chars_extracted = %.4f\n\n",
            (tot_time * 1000) / tot_ext);
    fprintf(stderr, "(Load_time+Extract_time)/Num_chars_extracted = %.4f\n\n",
            ((Load_time+tot_time) * 1000) / tot_ext);
}

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
    fprintf(stderr, "Usage:  %s <index> <type> [MAX_PAT] [SKIP]\n", progname);
    fprintf(stderr, "\n\t<type>   denotes the type of queries:\n");
    fprintf(stderr, "\t         %c counting queries;\n", COUNT);
    fprintf(stderr, "\t         %c locating queries;\n", LOCATE);
    fprintf(stderr, "\t         %c extracting queries.\n\n", EXTRACT);
    fprintf(stderr, "\n\t[MAX_PAT=100000000] Maximal MAX_PAT patterns are processed for count or locate.\n");
    fprintf(stderr, "\n\t[SKIP=0] The first SKIP patterns are skipped for count or locate.\n");
}
