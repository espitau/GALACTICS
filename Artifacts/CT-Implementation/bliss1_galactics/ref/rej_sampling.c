#include "randombytes.h"
#include "bliss.h"
#include "rej_sampling.h"
#if 0
#include <assert.h>
#else
#define assert(p) do{}while(0)
#endif

/* 
 * Computes the function:
 *
 *   f(x) = 2^41 * exp(-x/(2 sigma^2 * 256))
 *
 * over the interval [0, 2sigma^2 log(2) * 256] with over 39 bits of precision.
 */
uint64_t exp_scaled(uint64_t x)
{
    assert(x<=64082*256);

    uint64_t r;
    r = (( 809438661408LL    )*x) >> 28;
    r = (( 869506949331LL - r)*x) >> 28;
    r = (( 640044208952LL - r)*x) >> 27;
    r = (( 793458686015LL - r)*x) >> 27;
    r = (( 839743192604LL - r)*x) >> 27;
    r = (( 740389683060LL - r)*x) >> 26;
    r = ((1044449863563LL - r)*x) >> 27;
    r = (( 552517269260LL - r)*x) >> 25;
    r = (( 779422325990LL - r)*x) >> 23;
    r =  (2199023255552LL - r);

    return r;
}

/* 
 * Computes the function:
 *
 *   f(x) = 2^41 * exp(-x/(2 sigma^2))
 *
 * over the interval [0, 2sigma^2 log(2)] with over 39 bits of precision.
 */
uint64_t exp_simple(uint64_t x)
{
    return exp_scaled(x << 8);
}

/*
 * Computes the function:
 *
 *   f(x) = exp(x/sigma^2)
 *
 * over the interval [-B_2 sigma, B_2 sigma] with over 38 bits of
 * precision. Returns a 41-bit value r and fills the exponent variable
 * in such a way that (within >= 38 bits of precision)
 *
 *   f(x) = r / 2^exponent
 */
static uint64_t exp_extended(int64_t x, int *exponent)
{
    assert(x>-2787567 && x<2787567);

    uint64_t absx, y, k, r, b;
    int signx = (x >= 0);

    absx = CMUX(x, -x, (int64_t)signx) << 9;

    /* 
     * Constant-time Euclidean division by 
     * round(2sigma^2 log(2) * 256) = 16404853.
     * Note that the quotient fits within 7 bits.
     */
    k = 0;
    y = absx >> 7;
    for(int i = 6; i >= 0; i--) {
	y <<= 1;
	y  |= (absx >> i) & 1;
	b   = LSBMASK(y >= 16404853);
	k  |= (1UL << i) & b;
	y  -= 16404853 & b;
    }

    y = CMUX(16404853 - y, y, (uint64_t)signx);
    k+= signx;
    
    *exponent = CMUX(41 - k, 41 + k, signx);

    /*
     * x = 16404853 * (+-k) - y, hence:
     *
     * f(x) = exp(x / (2sigma^2 * 256)) 
     *      = exp(-y / (2sigma^2 * 256)) * exp(16404853/(2sigma^2*256))^(+-k)
     *
     * Now exp(16404853/(2sigma^2*256)) is very close to 2 by
     * construction; in fact, it is 2*(1 + 17933/2^43 + o(2^-46)). Since
     * k is at most 7 bits, the first-order expansion below ensures 39
     * bits of relative precision.
     */
    r = exp_scaled(y);
    k = (17933*k*r) >> 43;
    r = CMUX(r + k, r - k, (uint64_t)signx);

    return r;
}

/*
 * Computes the function:
 *
 *   f(x) = cosh(x/sigma^2)
 *
 * over the interval [-B_2 sigma, B_2 sigma] with over 38 bits of
 * precision. Returns a 41-bit value r and fills the exponent variable
 * in such a way that (within >= 38 bits of precision)
 *
 *   f(x) = r / 2^exponent
 */
static uint64_t cosh_simple(int64_t x, int *exponent)
{
    uint64_t r, r2;
    int64_t y = CABS(x);
    int exp2;

    r  = exp_extended( y, exponent);
    r2 = exp_extended(-y, &exp2);

    exp2 -= *exponent;
    r2 = CMUX(r2 >> exp2, 0, exp2<=40); 
    r += r2;

    *exponent += 1;

    return r;
}

unsigned char rejsamp_exp(uint64_t x)
{
    unsigned char rand[5];
    uint64_t p, u;

    p = exp_simple(x) >> 1;
    randombytes(rand, 5);
    u = rand[0] + ((uint64_t)rand[1] <<  8) + ((uint64_t)rand[2] << 16) + 
	          ((uint64_t)rand[3] << 24) + ((uint64_t)rand[4] << 32);

    return u >= p;
}

unsigned char rejsamp_cosh(int64_t x)
{
    unsigned char rand[5], r;
    uint64_t p, u0, u1, pu0, pu1;
    int exponent;

    p = cosh_simple(x, &exponent);
    randombytes(rand, 5);

    u0 = rand[0] + ((uint64_t)rand[1] << 8) +
	           ((uint64_t)(rand[2] & 0x0F) << 16);
    u1 = (rand[2] >> 4) + 
	           ((uint64_t)rand[3] << 4) + 
		   ((uint64_t)rand[4] << 12);
    u0 *= p;
    u1 *= p;
    pu0 = (u0 & 0xFFFFF) + ((u1 & 0x3FF) << 20);
    pu1 = (u0 >> 40) + (u1 >> 20) + (pu0 >> 40);
    pu0&= 0xFFFFF;

    /*
    if(exponent >= 0) {
	r = (pu1 >= (1ULL << exponent));
    }
    else if(exponent >= -40) {
	r = (pu0 >= (1ULL << (exponent+40))) || (pu1 != 0);
    }
    else {
	r = (pu1 != 0) || (pu0 != 0);
    }
    */

    unsigned char a, b, c, d;
    a = (pu1 >= (1ULL << exponent));
    b = (pu0 >= (1ULL << (exponent+40)));
    c = (pu1 != 0);
    d = (pu0 != 0);

    r = CMUX(CMUX(a, b||c, exponent >= 0),
	     c||d, exponent >= -40);

    return r;
}

