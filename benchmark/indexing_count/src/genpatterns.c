
// Extracts random patterns from a file

#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

static int Seed;
#define ACMa 16807
#define ACMm 2147483647
#define ACMq 127773
#define ACMr 2836
#define hi (Seed / ACMq)
#define lo (Seed % ACMq)

static int fst = 1;

/*
 * returns a random integer in 0..top-1
 */

int
aleat(top)
{
    long test;
    struct timeval t;
    if (fst) {
        gettimeofday(&t, NULL);
        Seed = t.tv_sec * t.tv_usec;
        fst = 0;
    }
    {
        Seed = ((test =
                     ACMa * lo - ACMr * hi) > 0) ? test : test + ACMm;
        return ((double) Seed) * top / ACMm;
    }
}

void parse_forbid(unsigned char* forbid, unsigned char** forbide)
{

    int len, i, j;
    len = strlen(forbid);

    *forbide = (unsigned char*) malloc((len+1)*sizeof(unsigned char));
    if (*forbide == NULL) {
        fprintf(stderr, "Error: cannot allocate %i bytes\n", len+1);
        fprintf(stderr, "errno = %i\n", errno);
        exit(1);
    }

    for (i = 0, j = 0; i < len; i++) {
        if (forbid[i] != '\\') {
            if (forbid[i] != '\n')
                (*forbide)[j++] = forbid[i];
        } else {
            i++;
            if (i == len) {
                forbid[i-1] = '\0';
                (*forbide)[j] = '\0';
                fprintf(stderr, "Not correct forbidden string: only one \\\n");
                return;
            }
            switch (forbid[i]) {
                case'n': (*forbide)[j++] = '\n'; break;
                case'\\': (*forbide)[j++] = '\\'; break;
                case'b': (*forbide)[j++] = '\b'; break;
                case'e': (*forbide)[j++] = '\e'; break;
                case'f': (*forbide)[j++] = '\f'; break;
                case'r': (*forbide)[j++] = '\r'; break;
                case't': (*forbide)[j++] = '\t'; break;
                case'v': (*forbide)[j++] = '\v'; break;
                case'a': (*forbide)[j++] = '\a'; break;
                case'c':
                    if (i+3 >= len) {
                        forbid[i-1] = '\0';
                        (*forbide)[j] = '\0';
                        fprintf(stderr, "Not correct forbidden string: 3 digits after \\c\n");
                        return;
                    }
                    (*forbide)[j++] = (forbid[i+1]-48)*100 +
                                      (forbid[i+2]-48)*10 + (forbid[i+3]-48);
                    i+=3;
                    break;
                default:
                    fprintf(stdout, "Unknown escape sequence '\\%c'in forbidden string\n", forbid[i]);
                    break;
            }
        }
    }
    (*forbide)[j] = '\0';
}

main(int argc, char** argv)
{
    int n, m, J, t;
    struct stat sdata;
    FILE* ifile, *ofile;
    unsigned char* buff;
    unsigned char* forbid, *forbide = NULL;

    if (argc < 5) {
        fprintf(stderr,
                "Usage: genpatterns <file> <length> <number> <patterns file> <forbidden>\n"
                "  randomly extracts <number> substrings of length <length> from <file>,\n"
                "  avoiding substrings containing characters in <forbidden>.\n"
                "  The output file, <patterns file> has a first line of the form:\n"
                "    # number=<number> length=<length> file=<file> forbidden=<forbidden>\n"
                "  and then the <number> patterns come successively without any separator.\n"
                "  <forbidden> uses \\n, \\t, etc. for nonprintable chracters or \\cC\n"
                "  where C is the ASCII code of the character written using 3 digits.\n\n");
        exit(1);
    }

    if (stat(argv[1], &sdata) != 0) {
        fprintf(stderr, "Error: cannot stat file %s\n", argv[1]);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }
    n = sdata.st_size;

    m = atoi(argv[2]);
    if ((m <= 0) || (m > n)) {
        fprintf(stderr,
                "Error: length must be >= 1 and <= file length"
                " (%i)\n", n);
        exit(1);
    }

    J = atoi(argv[3]);
    if (J < 1) {
        fprintf(stderr, "Error: number of patterns must be >= 1\n");
        exit(1);
    }

    if (argc > 5) {
        forbid = argv[5];
        parse_forbid(forbid, &forbide);
    } else
        forbid = NULL;

    ifile = fopen(argv[1], "r");
    if (ifile == NULL) {
        fprintf(stderr, "Error: cannot open file %s for reading\n", argv[1]);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }

    buff = (unsigned char*) malloc(n);
    if (buff == NULL) {
        fprintf(stderr, "Error: cannot allocate %i bytes\n", n);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }

    if (fread(buff, n, 1, ifile) != 1) {
        fprintf(stderr, "Error: cannot read file %s\n", argv[1]);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }
    fclose(ifile);

    ofile = fopen(argv[4], "w");
    if (ofile == NULL) {
        fprintf(stderr, "Error: cannot open file %s for writing\n",
                argv[4]);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }

    if (fprintf(ofile, "# number=%i length=%i file=%s forbidden=%s\n",
                J, m, argv[1],
                forbid == NULL ? "" : (char*) forbid) <= 0) {
        fprintf(stderr, "Error: cannot write file %s\n", argv[4]);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }

    for (t = 0; t < J; t++) {
        int j, l;
        if (!forbide)
            j = aleat(n - m + 1);
        else {
            do {
                j = aleat(n - m + 1);
                for (l = 0; l < m; l++)
                    if (strchr(forbide, buff[j + l]))
                        break;
            } while (l < m);
        }
        for (l = 0; l < m; l++)
            if (putc(buff[j + l], ofile) != buff[j + l]) {
                fprintf(stderr,
                        "Error: cannot write file %s\n",
                        argv[4]);
                fprintf(stderr, " errno = %i\n", errno);
                exit(1);
            }
    }

    if (fclose(ofile) != 0) {
        fprintf(stderr, "Error: cannot write file %s\n", argv[4]);
        fprintf(stderr, " errno = %i\n", errno);
        exit(1);
    }

    fprintf(stderr, "File %s successfully generated\n", argv[4]);
    free(forbide);
    exit(0);
}
