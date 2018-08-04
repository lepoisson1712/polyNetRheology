#ifndef RNGMIT_H_
#define RNGMIT_H_

/* for 64-bit integers: #define rng_conv 5.421010862427522e-20 */
/* for 32-bit integers: #define rng_conv 2.3283064365387e-10   */
#if defined(__LP64__) || defined(_LP64)
#define BIT64   1
#endif
/* May 21 2010: try to do this automatically at compile time */

#ifdef BIT64	
#define rng_conv 5.421010862427522e-20
#else
#define rng_conv 2.3283064365387e-10
#endif

#define rngmit (rng_conv*(rng_ia[rng_p=rng_mod[rng_p]] += rng_ia[rng_pp=rng_mod[rng_pp]]))
#define rngmitint (rng_ia[rng_p=rng_mod[rng_p]] += rng_ia[rng_pp=rng_mod[rng_pp]])

extern unsigned long int rng_ia[55];
extern int rng_p,rng_pp;
extern int rng_mod[55];

void rngseed(unsigned long int s);
void rnginit (unsigned long int ia[], int pp, int p);
#endif // RNGMIT_H_
