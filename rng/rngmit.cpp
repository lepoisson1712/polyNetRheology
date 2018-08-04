#include <stdio.h>
#include "rngmit.hpp"

/* Implementation of the additive random number generator of Mitchell and
 * Moore.  The generator itself is implemented as a macro in the accompanying
 * header file.
 *
 * Mark Newman  8 OCT 96
 */

/* Constants */

#define a 2416
#define c 374441
#define m 1771875
#define conv 2423.96743336861

/* Globals */

static unsigned long int i;
unsigned long int rng_ia[55];
int rng_p,rng_pp;
int rng_mod[55];

/* Function to seed all the random number generators */
void rngseed(unsigned long int s)
{
  int n;

  /* First seed the linear congruential generator */

#ifdef BIT64
  printf("# 64-BIT VERSION: rng_conv=%e\n",rng_conv);
#else
  printf("# 32-BIT VERSION: rng_conv=%e\n",rng_conv);
#endif

  i = s;

  /* Use that to seed the additive generator.  Also setup the mod array */

  for (n=0; n<55; n++) {
    rng_ia[n] = (long unsigned int) (conv*(i=(a*i+c)%m));
    rng_mod[n] = n-1;
  }
  rng_mod[0] = 54;

  rng_p = 0;
  rng_pp = 24;

  /* Run off ten thousand random numbers, just to get things started */

  for (n=0; n<10000; n++) rng_ia[rng_p=rng_mod[rng_p]] += rng_ia[rng_pp=rng_mod[rng_pp]];
}

void rnginit (unsigned long int ia[], int pp, int p)
{
	  for (int n=0; n<55; n++) {
	    rng_ia[n] = ia[n];
	    rng_mod[n] = n-1;
	  }
	  rng_mod[0] = 54;

	  rng_p = p;
	  rng_pp = pp;
}
