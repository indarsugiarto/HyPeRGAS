#include <sark.h>

/* Indar: since we want to generate "float" number */
#include <stdfix.h>
#include <math.h>

#define REAL            accum
#define REAL_CONST(x)   x##k
#define RAND_MAX_UINT	0xFFFFFFFF
#define RAND_MAX_UINT16	0xFFFF

/* Period parameters */  
#define N			624
#define M			397
#define MATRIX_A	0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK	0x80000000UL /* most significant w-r bits */
#define LOWER_MASK	0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
uint genrand_int32(void)
{
    uint rn;
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    rn = y;
    return rn;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* genrand_fixp(min, max) generates random values within [min, max]
 *
 * Source: http://stackoverflow.com/questions/7978759/generate-float-random-values-also-negative
 * In general to generate random numbers from an arbitrary distribution you'd first generate uniform random numbers and then pass them to the inverse of the cumulative distribution function.

Assume for example that you want random numbers with uniform distribution on the interval [-10.0, 10.0] and all you've got is random numbers from [0.0, 1.0]. Cumulative distribution function of the uniform distribution on [-10.0, 10.0] is:

cdf(x) = 0.05 * x + 0.5   for x in [-10.0, 10.0]
This expresses the probability that a random number generated is smaller than x. The inverse is

icdf(y) = 20.0 * y - 10.0   for y in [0.0, 1.0]
(You can obtain this easily on paper by switching the x and y axis).

Hence to obtain random numbers uniformly distributed on [-10.0, 10.0] you can use the following code:

#include <stdlib.h>

// Returns uniformly distributed random numbers from [0.0, 1.0].
double uniform0to1Random() {
	double r = random();
	return r / ((double)RAND_MAX + 1);
}

// Returns uniformly distributed random numbers from [-10.0, 10.0].
double myRandom() {
  return 20.0 * uniform0to1Random() - 10.0;
}
In fact, you don't need uniform0to1Random() since there are already a lot of good uniform random numbers generators from [0.0, 1.0] (e.g. in the boost library).

You can use the method to generate random numbers with nearly any probability distribution you want by sampling the inverse cumulative distribution as shown above.

 *
 * */
REAL genrand_fixp(REAL minVal, REAL maxVal, uint seed)
{
	//init_genrand(seed);
	ushort rr = genrand_int32() >> 16;
	REAL upper = REAL_CONST(0.0);
	REAL ratio = (maxVal-minVal)/maxVal;
	float f = (float)rr / ((float)RAND_MAX_UINT16);
	REAL R = f;
	//REAL R = (REAL)rr / ((REAL)RAND_MAX_UINT16 + upper);
	REAL result = ratio * maxVal * R + minVal;
	return result;
}


