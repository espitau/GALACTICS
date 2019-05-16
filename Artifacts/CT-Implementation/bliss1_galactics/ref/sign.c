#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "api.h"
#include "crypto_sign.h"
#include "randombytes.h"
#include "fips202.h"
#include "poly.h"
#include "reduce.h"
#include "gaussian_sampling.h"
#include "rej_sampling.h"
#include "packing.h"

#define generate_c_consttime generate_c

static void compute_u(poly_t *u, const polysmall_t *y1, const polysmall_t *y2,
	const poly_t *a_fft);

static void generate_c_vartime(hpol_t *c, poly_t *c_expand,
               const unsigned char mu[CRHBYTES],
               const poly_t *u)
{
  unsigned int i, b, pos;
  unsigned char inbuf[PACK16_BYTES + CRHBYTES];
  unsigned char outbuf[SHAKE128_RATE];
  uint64_t state[25];
  
  poly_pack16(inbuf, u);
  memcpy(inbuf + PACK16_BYTES, mu, CRHBYTES);

  shake128_absorb(state, inbuf, PACK16_BYTES + CRHBYTES);
  shake128_squeezeblocks(outbuf, 1, state);

  pos = 0;

  for(i = 0; i < N; ++i)
    c_expand->coeffs[i] = 0;

  for(i = N-kappa; i < N; i++) {
    do {
      if(pos >= SHAKE128_RATE) {
        shake128_squeezeblocks(outbuf, 1, state);
        pos = 0;
      }

      b = outbuf[pos++];
      b = (b&1) << 8;
      b|= outbuf[pos++];
    } while(b > i);

    c_expand->coeffs[i] = c_expand->coeffs[b];
    c_expand->coeffs[b] = 1;
  }
  
  pos = 0;
  for(i = 0; i < kappa; i++)
      c->indices[i] = 0;

  for(i = 0; i < N; i++) {
      if(c_expand->coeffs[i])
	  c->indices[pos++] = i;
  }

}

static void generate_c_consttime(hpol_t *c, poly_t *c_expand,
               const unsigned char mu[CRHBYTES],
               const poly_t *u)
{
    unsigned int i, j, b, pos, cmp, t;
    unsigned char inbuf[PACK16_BYTES + CRHBYTES];
    unsigned char outbuf[SHAKE128_RATE];
    uint64_t state[25];
  
    poly_pack16(inbuf, u);
    memcpy(inbuf + PACK16_BYTES, mu, CRHBYTES);

    shake128_absorb(state, inbuf, PACK16_BYTES + CRHBYTES);
    shake128_squeezeblocks(outbuf, 1, state);

    pos = 0;

    for(j=0; j < kappa; j++)
	c->indices[j] = 0;

    for(i = N-kappa; i < N; i++) {
	do {
	if(pos >= SHAKE128_RATE) {
	    shake128_squeezeblocks(outbuf, 1, state);
	    pos = 0;
	}

	b = outbuf[pos++];
	b = (b&1) << 8;
	b|= outbuf[pos++];
	} while(b > i);

	for(j=0; j + N < i + kappa; j++) {
	    cmp = -(c->indices[j] > b);
	    t   = c->indices[j] ^ b;
	    b  ^= t & cmp;
	    b  ^= (-!t) & (b^i);

	    c->indices[j] ^= t & cmp;
	}
	c->indices[j] = b;
    }
    /*
    printf("c = ");
	for(int j=0; j<kappa; j++)
	    printf("%" PRIu16 "%c", c->indices[j], (j<kappa-1)?' ':'\n');
    */

    for(i=0; i<N; i++) {
	c_expand->coeffs[i] = 0;
	for(j=0; j<kappa; j++) 
	    c_expand->coeffs[i] ^= (c->indices[j] == i);
    }

}


static int32_t normNkS(const polysmall_t *restrict s1,
	const polysmall_t *restrict s2);

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    polysmall_t s1, s2;
    poly_t invs1, a_fft;
    poly_t s1_fft, s2_fft;

    do {
	poly_sample_sparse(&s1);
	poly_sample_sparse(&s2);

	for(int i=0; i<N; i++)
	    s2.coeffs[i] = 2*s2.coeffs[i];
	s2.coeffs[0] += 1;

	if(normNkS(&s1, &s2) > boundNkS)
	    continue;

	poly_unsmall(&s1_fft, &s1);
	poly_ntt(&s1_fft);
	poly_freeze(&s1_fft);

	for(int i=0; i<N; i++) {
	    if(s1_fft.coeffs[i] == 0)
		continue;
	    invs1.coeffs[i] = inverse_modq(s1_fft.coeffs[i]);
	}

	poly_unsmall(&s2_fft, &s2);
	poly_ntt(&s2_fft);
	poly_freeze(&s2_fft);
	
	for(int i=0; i<N; i++) {
	    a_fft.coeffs[i] = freeze(s2_fft.coeffs[i] * invs1.coeffs[i]);
	}

	poly_t u;
	compute_u(&u, &s1, &s2, &a_fft);
	if(u.coeffs[0] != Q)
	    continue;
	for(int i=1; i<N; i++)
	    if(u.coeffs[i] != 0)
		continue;
	break;
    } while(1);

    pack_pk(pk, &a_fft);
#if MULTC_FFT
    pack_sk_fft(sk, pk, &s1_fft, &s2_fft);
#else
    pack_sk(sk, pk, &s1, &s2);
#endif

    return 0;
}

static inline int64_t innerprod(const polysmall_t *x, const polysmall_t *y) {
    int64_t r = 0, xi, yi;

    for(int i=0; i<N; i++) {
	xi = x->coeffs[i];
	yi = y->coeffs[i];
	r += xi*yi;
    }

    return r;
}

static inline int64_t norm2(const polysmall_t *x) {
    return innerprod(x, x);
}

static inline int32_t norminf(const polysmall_t *x) {
    int32_t r = 0, xi;

    for(int i=0; i<N; i++) {
	xi = x->coeffs[i]; 
	r  = CMUX( xi, r,  xi > r);
	r  = CMUX(-xi, r, -xi > r);
    }

    return r;
}

static inline int boundcheck_nocomp(const polysmall_t *z1, const polysmall_t *z2) {
    int64_t norm2z;
    int32_t norminfz1, norminfz2;
    int r;

    norm2z    = norm2(z1) + norm2(z2);
    norminfz1 = norminf(z1);
    norminfz2 = norminf(z2);
    
    r = (norm2z >= B_VERIFY);
    r = (norminfz1 >= BINFTY) || r;
    r = (norminfz2 >= BINFTY) || r;

    return r;
}

static inline int boundcheck(const polysmall_t *z1, const polysmall_t *z2dag) {
    int64_t norm2z;
    int32_t norminfz1, norminfz2;
    int r;

    norm2z    = norm2(z1) + (norm2(z2dag) << (2*dropped_bits));
    norminfz1 = norminf(z1);
    norminfz2 = norminf(z2dag) << dropped_bits;
    
    r = (norm2z >= B_VERIFY);
    r = (norminfz1 >= BINFTY) || r;
    r = (norminfz2 >= BINFTY) || r;

    return r;
}

static void compute_u(poly_t *u, const polysmall_t *y1, const polysmall_t *y2,
	const poly_t *a_fft)
{
    poly_unsmall(u, y1);
    
    poly_ntt(u);
    poly_pointwise_invmontgomery(u, u, a_fft);
    poly_invntt_montgomery(u);

    for(int i=0; i<N; i++) {
	u->coeffs[i] = reduce2q(2*u->coeffs[i]*oneQTwo + 2*Q + y2->coeffs[i]);
    }
}

static void round_u(poly_t *v, const poly_t *u)
{
    for(int i=0; i<N; i++) {
	/*
	uint32_t vi = u->coeffs[i];
	int sign = (vi >= Q);

	vi = CMUX(vi, 2*Q-vi, sign);
	vi = 2*vi + (1<<dropped_bits);
	vi = reducep(vi >> (dropped_bits + 1));

	v->coeffs[i] = CMUX(vi, modp-vi, sign);
	*/
	uint32_t vi = u->coeffs[i];
	vi = 2*vi + (1<<dropped_bits);
	vi = reducep(vi >> (dropped_bits + 1));
	v->coeffs[i] = vi;
    }
}

static void mult_by_c(polysmall_t *sc, const polysmall_t *s, const hpol_t *c)
{
#if 0
    int32_t sci;
    int ind;
    for(int i=0; i<N; i++) {
	sci = 0;

	/* Note that this code leaks c via access patterns,
	 * but this is ok! */
	for(int j=0; j<kappa; j++) {
	    ind = i-c->indices[j];
	    sci += (ind>=0) ? s->coeffs[ind] : -s->coeffs[N+ind];
	}
	sc->coeffs[i] = sci;
    }
#elif 1
    polysmall_t rot;
    for(int i=0; i<N; i++) {
	sc->coeffs[i] = 0;
    }
    for(int i=0; i<kappa; i++) {
	poly_rot(&rot, s, c->indices[i]);
	poly_add(sc, &rot);
    }

#else
    poly_t sbig, cbig;
    poly_unsmall(&sbig, s);

    for(int i=0; i<N; i++)
	cbig.coeffs[i] = 0;
    for(int i=0; i<kappa; i++)
	cbig.coeffs[c->indices[i]] = 1;
    poly_ntt(&sbig);
    poly_ntt(&cbig);
    poly_pointwise_invmontgomery(&sbig, &sbig, &cbig);
    poly_invntt_montgomery(&sbig);
    poly_freeze(&sbig);
    poly_small(sc, &sbig);
#endif
}

static void compress_z2(polysmall_t *z2dag, const polysmall_t *z2, const poly_t *u, 
	const poly_t *ured)
{
    poly_t umx;
    int32_t z2i;

    for(int i=0; i<N; i++) {
	umx.coeffs[i] = reduce2q(u->coeffs[i] + 2*Q - z2->coeffs[i]);
    }

    round_u(&umx, &umx);

    for(int i=0; i<N; i++) {
	z2i = reducep(ured->coeffs[i] + 2*modp - umx.coeffs[i]);
	z2dag->coeffs[i] = CMUX(z2i, z2i-modp, 2*z2i < modp);
    }
}

int crypto_sign(unsigned char *sm,
                unsigned long long *smlen,
                const unsigned char *m,
                unsigned long long mlen,
                const unsigned char *sk)
{
    polysmall_t z1, z2, s1c, s2c;
    poly_t a_fft, u, c_fft;
   
#if MULTC_FFT
    poly_t s1_fft, s2_fft, sc_fft;
#else
    polysmall_t s1, s2;
#endif
    hpol_t c;
#if COMPRESS
    polysmall_t z2dag;
    poly_t ured;
#endif
    unsigned char reject = 1, rejexp, rejcosh, rejbound;

    unsigned char hash[CRHBYTES], bit;

#if MULTC_FFT
    unpack_sk_fft(&s1_fft, &s2_fft, &a_fft, sk);
#else
    unpack_sk(&s1, &s2, &a_fft, sk);
#endif
    shake128(hash, CRHBYTES, m, mlen);

    *smlen = CRYPTO_BYTES + mlen;
    memmove(sm + CRYPTO_BYTES, m, mlen);

    while(reject) {
	sample_gaussian_poly(&z1, &z2);
	compute_u(&u, &z1, &z2, &a_fft);
#if COMPRESS
	round_u(&ured, &u);
	generate_c(&c, &c_fft, hash, &ured);
#else
	generate_c(&c, &c_fft, hash, &u);
#endif
#if MULTC_FFT
	poly_ntt(&c_fft);
	poly_pointwise_invmontgomery(&sc_fft, &s1_fft, &c_fft);
	poly_invntt_montgomery(&sc_fft);
	poly_small(&s1c, &sc_fft);
	poly_pointwise_invmontgomery(&sc_fft, &s2_fft, &c_fft);
	poly_invntt_montgomery(&sc_fft);
	poly_small(&s2c, &sc_fft);
#else	
	mult_by_c(&s1c, &s1, &c);
	mult_by_c(&s2c, &s2, &c);
#endif
	randombytes(&bit, 1);
	poly_sign_flip(&s1c, bit);
	poly_sign_flip(&s2c, bit);

	poly_add(&z1, &s1c);
	poly_add(&z2, &s2c);

#if COMPRESS
	compress_z2(&z2dag, &z2, &u, &ured);
#endif

	rejexp   = rejsamp_exp(norm2(&s1c) + norm2(&s2c));
	rejcosh  = rejsamp_cosh(innerprod(&s1c,&z1) + innerprod(&s2c,&z2));
#if COMPRESS
	rejbound = boundcheck(&z1, &z2dag);
#else
	rejbound = boundcheck_nocomp(&z1, &z2);
#endif

	reject   = rejcosh  || rejexp;
	reject   = rejbound || reject;
    }

#if COMPRESS
    pack_sig(sm, &z1, &z2dag, &c);
#else
    pack_sig_nocomp(sm, &z1, &z2, &c);
#endif
    return 0;
}

int crypto_sign_open(unsigned char *m,
                     unsigned long long *mlen,
                     const unsigned char *sm,
                     unsigned long long smlen,
                     const unsigned char *pk)
{
    unsigned char hash[CRHBYTES];
  
    poly_t a_fft, u;
    polysmall_t z1;
#if COMPRESS
    polysmall_t z2dag, cpol;
#else
    polysmall_t z2;
#endif
    
    hpol_t c, csig;
  
    if(smlen < CRYPTO_BYTES)
      goto badsig;
  
    *mlen = smlen - CRYPTO_BYTES;
  
    unpack_pk(&a_fft, pk);

#if COMPRESS
    if(unpack_sig(&z1, &z2dag, &csig, sm)) {
	goto badsig;
    }
    if(boundcheck(&z1, &z2dag)) {
	goto badsig;
    }
#else
    if(unpack_sig_nocomp(&z1, &z2, &csig, sm)) {
	goto badsig;
    }
    if(boundcheck_nocomp(&z1, &z2)) {
	goto badsig;
    }
#endif 

    shake128(hash, CRHBYTES, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);

#if COMPRESS
    for(int i=0; i<N; i++)
        cpol.coeffs[i] = 0;
    for(int i=0; i<kappa; i++)
        cpol.coeffs[csig.indices[i]] = Q*oneQTwo;
  
    compute_u(&u, &z1, &cpol, &a_fft);
    round_u(&u, &u);
  
    for(int i=0; i<N; i++) {
        u.coeffs[i] += modp + z2dag.coeffs[i];
        u.coeffs[i] %= modp;
    }
#else
    compute_u(&u, &z1, &z2, &a_fft);
    for(int i=0; i<kappa; i++) {
        u.coeffs[csig.indices[i]] += Q;
        u.coeffs[csig.indices[i]] %= 2*Q;
    }
#endif

    poly_t c_expand;
    generate_c_vartime(&c, &c_expand, hash, &u);

    for(int i=0; i<kappa; i++) {
	if(c.indices[i] != csig.indices[i]) {
	    goto badsig;
	}
    }

    memmove(m, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
    return 0;
  
badsig:
    /* Signature verification failed */
    *mlen = (unsigned long long) -1;
    memset(m, 0, smlen - CRYPTO_BYTES);

    return -1;
}

int cmpdec(const void *a, const void *b) {
    const int32_t *A = a, *B = b;
    return (*A < *B) - (*A > *B);
}

static int32_t normNkS(const polysmall_t *restrict s1,
	const polysmall_t *restrict s2)
{
    polysmall_t gram, rot;
    int32_t max[N], r = 0;

    for(int i=0; i<N; i++) {
	poly_rot(&rot, s1, i);
	gram.coeffs[i] = innerprod(&rot, s1);
	poly_rot(&rot, s2, i);
	gram.coeffs[i]+= innerprod(&rot, s2);
    }

    for(int i=0; i<N; i++) {
	poly_rot(&rot, &gram, i);
	qsort(rot.coeffs, N, sizeof(int32_t), cmpdec);
	max[i] = 0;
	for(int j=0; j<kappa; j++)
	    max[i] += rot.coeffs[j];
    }

    qsort(max, N, sizeof(int32_t), cmpdec);
    for(int j=0; j<kappa; j++)
	r += max[j];

    return r;

}
