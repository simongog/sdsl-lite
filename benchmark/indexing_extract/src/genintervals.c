
// Chooses random intervals from a file

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

static int fst=1;

/* returns a random integer in 0..top-1 */

int aleat(top)

{
    long test;
    struct timeval t;
    if (fst) {
        gettimeofday(&t,NULL);
        Seed = t.tv_sec*t.tv_usec;
        fst=0;
    }
    {
        Seed = ((test = ACMa * lo - ACMr * hi) > 0) ? test : test + ACMm;
        return ((double)Seed) * top/ACMm;
    }
}

main(int argc, char** argv)

{
    int n,m,J,t;
    struct stat sdata;
    FILE* ifile,*ofile;
    unsigned char* buff;

    if (argc < 5) {
        fprintf(stderr,
                "Usage: %s <file> <length> <number> <intervals file>\n"
                "  randomly extracts <number> intervals of length <length> from <file>.\n"
                "  The output file, <intervals file> has a first line of the form:\n"
                "    # number=<number> length=<length> file=<file>\n"
                "  and then <number> lines of the form <from>,<to>.\n",argv[0]
               );
        exit(1);
    }

    if (stat(argv[1],&sdata) != 0) {
        fprintf(stderr,"Error: cannot stat file %s\n",argv[1]);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }
    n = sdata.st_size;

    m = atoi(argv[2]);
    if ((m <= 0) || (m > n)) {
        fprintf(stderr,"Error: length must be >= 1 and <= file length"
                " (%i)\n",n);
        exit(1);
    }

    J = atoi(argv[3]);
    if (J < 1) {
        fprintf(stderr,"Error: number of intervals must be >= 1\n");
        exit(1);
    }

    ifile = fopen(argv[1],"r");
    if (ifile == NULL) {
        fprintf(stderr,"Error: cannot open file %s for reading\n",argv[1]);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }

    buff = (unsigned char*) malloc(n);
    if (buff == NULL) {
        fprintf(stderr,"Error: cannot allocate %i bytes\n",n);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }

    if (fread(buff,n,1,ifile) != 1) {
        fprintf(stderr,"Error: cannot read file %s\n",argv[1]);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }
    fclose(ifile);

    ofile = fopen(argv[4],"w");
    if (ofile == NULL) {
        fprintf(stderr,"Error: cannot open file %s for writing\n",argv[4]);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }

    if (fprintf(ofile,"# number=%i length=%i file=%s\n", J,m,argv[1]) <= 0) {
        fprintf(stderr,"Error: cannot write file %s\n",argv[4]);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }

    for (t=0; t<J; t++) {
        int j,l;
        j = aleat(n-m+1);
        if (fprintf(ofile,"%i,%i\n",j,j+m-1) <= 0) {
            fprintf(stderr,"Error: cannot write file %s\n",argv[4]);
            fprintf(stderr," errno = %i\n",errno);
            exit(1);
        }
    }

    if (fclose(ofile) != 0) {
        fprintf(stderr,"Error: cannot write file %s\n",argv[4]);
        fprintf(stderr," errno = %i\n",errno);
        exit(1);
    }

    fprintf(stderr,"File %s successfully generated\n",argv[4]);
    exit(0);// SG: Replaced 1 by 0.
}


