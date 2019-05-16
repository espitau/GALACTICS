#include <inttypes.h>
#include <string.h>
#include "randombytes.h"
#include "bliss.h"
#include "gaussian_sampling.h"
#include "rej_sampling.h"

#define sample_gaussian_poly_bern sample_gaussian_poly 

static uint32_t sample_berncdt(uint64_t utop, uint64_t ubot)
{
    uint32_t x = 0;
    uint32_t b = 1;
    uint32_t r, s, t;
    t  = (utop >  11881272476311950404ULL);
    r  = (utop == 11881272476311950404ULL);
    s  = (ubot >   2232598800125794762ULL);
    b  = (r && s) || t;
    x += b;
    t  = (utop >  17729174351313943813ULL);
    r  = (utop == 17729174351313943813ULL);
    s  = (ubot >  17046599807202850264ULL);
    b  = (r && s) || t;
    x += b;
    t  = (utop >  18426461144592799266ULL);
    r  = (utop == 18426461144592799266ULL);
    s  = (ubot >   9031501729263515114ULL);
    b  = (r && s) || t;
    x += b;
    t  = (utop >  18446602887327906610ULL);
    r  = (utop == 18446602887327906610ULL);
    s  = (ubot >  11817852927693963396ULL);
    b  = (r && s) || t;
    x += b;
    t  = (utop >  18446743834670245612ULL);
    r  = (utop == 18446743834670245612ULL);
    s  = (ubot >   7306021935394802834ULL);
    b  = (r && s) || t;
    x += b;
    t  = (utop >  18446744073611412414ULL);
    r  = (utop == 18446744073611412414ULL);
    s  = (ubot >  17880792342251759005ULL);
    b  = (r && s) || t;
    x += b;
    t  = (utop >  18446744073709541852ULL);
    r  = (utop == 18446744073709541852ULL);
    s  = (ubot >  14689009182029885173ULL);
    b  = (r && s) || t;
    x += b;
    r  = (utop == 18446744073709551615ULL);
    s  = (ubot >  14106032229701791861ULL);
    b  = (r && s);
    x += b;
    s  = (ubot >  18446718728838181855ULL);
    b &= s;
    x += b;
    s  = (ubot >  18446744073673701140ULL);
    b &= s;
    x += b;
    
    return x;
}


static void sample_gaussian_poly_cdt(polysmall_t *y) 
{
    int64_t yi;

    uint64_t *us;
    unsigned char *cs;
    unsigned char buf[8*N*sizeof(uint64_t) + N];

    start_time("randombytes");
    randombytes(buf, sizeof buf);
    stop_time("randombytes");

    cs = buf;
    us = (uint64_t*)(buf+N);

    start_time("sample_mw_cdt");
    for(int i = 0; i < N; i++) {
	yi = sample_mw_cdt(us,*cs++);
	us+= 8;
	y->coeffs[i] = yi;
    }
    stop_time("sample_mw_cdt");
}

static inline uint32_t div16404853(uint32_t x) {
    uint32_t y, z;
    y = (97488647 * (uint64_t) x) >> 32;
    z = (x - y) >> 1;
    return (z + y) >> 23;
}

#define BERN_RANDMULT	1664

void sample_gaussian_poly_bern(polysmall_t *p1, polysmall_t *p2) 
{
    size_t i = 0, pos = 0;
    uint32_t x, y, t, k;
    int32_t z;

    uint64_t *us, utop, ubot, v, w;
    unsigned char *cs, c;
    unsigned char buf[(3*sizeof(uint64_t) + 2)*BERN_RANDMULT];
    uint32_t coeffs[2*N];

    randombytes(buf, sizeof buf);

    cs = buf;
    us = (uint64_t*)(buf+2*BERN_RANDMULT);

    while(i<2*N) {
	if(pos++ >= BERN_RANDMULT) {
	    randombytes(buf, sizeof buf);
	    pos = 0;
	    cs = buf;
	    us = (uint64_t*)(buf+2*BERN_RANDMULT);
	}

	utop = *us++;
	ubot = *us++;
	w    = *us++;
	w  >>= 1;

	x = sample_berncdt(utop, ubot) << 8;
	y = *cs++;
	z = x + y;

	y = y*(y + 2*x) << 8;  
	k = div16404853(y);
	t = y - 16404853*k;
	v = exp_scaled(t) << (22-k);
	
	c = *cs++;
	if((w > v) || (c & (z==0)))
	   continue;

	c >>= 1;
	coeffs[i++] = CFLIP(z,(int32_t)c);
    }
    for(int j=0; j<N; j++)
	p1->coeffs[j] = coeffs[j];
    for(int j=0; j<N; j++)
	p2->coeffs[j] = coeffs[N+j];
}
