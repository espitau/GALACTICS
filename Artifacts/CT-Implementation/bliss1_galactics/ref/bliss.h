#ifndef BLISS_H
#define BLISS_H

#include <inttypes.h>

/* Constant-time macros */
#define LSBMASK(c)	(-((c)&1))
#define	CMUX(x,y,c)	(((x)&(LSBMASK(c)))^((y)&(~LSBMASK(c))))
#define CFLIP(x,c)	CMUX(x,-(x),c)
#define CABS(x)		CFLIP(x,x>=0)

#define N		512
#define log2N		9
#define alpha_rejection	1
#define Q		12289
#define kappa		23
#define density		0.3
#define density2	0
#define	d1		154 // ceil(N*density)
#define	d2		0   // ceil(N*density2)
#define C		1.62
#define sigma		215
#define B_VERIFY	12872*12872
#define BINFTY		2100
#define dropped_bits	10
#define oneQTwo		18433 // 1/(q+2) mod 2q
#define modp		24
#define boundNkS	46478 // floor(C*C*5*(d1+4*d2)*kappa)

#define CRHBYTES	32



typedef struct {
    uint32_t	coeffs[N];
} poly_t;

typedef struct {
    int32_t	coeffs[N];
} polysmall_t;


typedef struct {
    uint16_t	indices[kappa];
} hpol_t;

#endif
