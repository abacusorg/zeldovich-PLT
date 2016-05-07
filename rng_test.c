/* The purpose of this test is to check the cross-platform
 * reproducibility of our random number generation.  This is
 * important if we want to generate IC files with the same
 * wave modes for a given random seed.
 * Usage: make run_rng_test
 */

#include <stdio.h>
#include <gsl/gsl_rng.h>

int main (void)
{
    const gsl_rng_type * T;
    gsl_rng * r;

    int Nstart = 1000000;
    int Nend = Nstart + 100;

    gsl_rng_env_setup();

    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    int i;
    for(i = 0; i < Nend; i++)
        if(i > Nstart)
            printf("%f\n", gsl_rng_uniform (r));

    gsl_rng_free (r);

    return 0;
}
