/* generate random sequences over a specified alphabet size, of
   a given length, and with a given number of repeats.
   each repeat commences with a different character, but then has
   all remaining characters identical

   Alistair Moffat, April 2013.
*/

#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char* argv[])
{

    int c;
    int sigma=64;
    int base=1000;
    int repeats=100;
    int i;

    if (argc==1) {
        fprintf(stderr, "Usage: %s sigma baselen repeats\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
    if (argc>1) {
        sigma = atoi(argv[1]);
        if (sigma > 90) {
            fprintf(stderr, "sigma of %d is too big\n",
                    sigma);
            exit(EXIT_FAILURE);
        }
    }
    if (argc>2) {
        base = atoi(argv[2]);
    }
    if (argc>3) {
        repeats = atoi(argv[3]);
    }
    fprintf(stderr, "sigma=%d, base=%d, repeats=%d\n",
            sigma, base, repeats);
    while (repeats) {
        srand(7716977);
        putchar(32+repeats);
        for (i=1; i<base; i++) {
            c = 32 + rand() % sigma;
            putchar(c);
        }
        repeats--;
    }
    return 0;
}
