#include <stdint.h>
#include <stddef.h>
#include "fips202.h"
#include "bliss.h"
#include "reduce.h"
#include "ntt.h"
#include "poly.h"
#include "randombytes.h"

void poly_sample_sparse(polysmall_t *a)
{
    unsigned char buf[512];
    int ind, num, sign;
    size_t pos = 0;

    for(int i=0; i<N; i++)
	a->coeffs[i] = 0;

    randombytes(buf, sizeof buf);
    for(num = 0; num < d1; ) {
	if(pos+1 >= sizeof buf) {
	    randombytes(buf, sizeof buf);
	    pos = 0;
	}
	ind  = buf[pos++];
	ind |= (buf[pos] & 1) << 8;
	sign = buf[pos++] >> 1;
	if(a->coeffs[ind] == 0) {
	    a->coeffs[ind] = CFLIP(1,sign);
	    num++;
	}
    }
    for(num = 0; num < d2; ) {
	if(pos+1 >= sizeof buf) {
	    randombytes(buf, sizeof buf);
	    pos = 0;
	}
	ind  = buf[pos++];
	ind |= (buf[pos] & 1) << 8;
	sign = buf[pos++] >> 1;
	if(a->coeffs[ind] == 0) {
	    a->coeffs[ind] = CFLIP(2,sign);
	    num++;
	}
    }
}

void poly_reduce(poly_t *a) {
  unsigned int i;

  for(i = 0; i < N; ++i)
    a->coeffs[i] = reduce32(a->coeffs[i]);

}

void poly_csubq(poly_t *a) {
  unsigned int i;

  for(i = 0; i < N; ++i)
    a->coeffs[i] = csubq(a->coeffs[i]);

}

void poly_freeze(poly_t *a) {
  unsigned int i;

  for(i = 0; i < N; ++i)
    a->coeffs[i] = freeze(a->coeffs[i]);

}

void poly_add(polysmall_t *restrict a, const polysmall_t *restrict b)  {
  unsigned int i;

  for(i = 0; i < N; ++i)
    a->coeffs[i] += + b->coeffs[i];

}

void poly_sign_flip(polysmall_t *a, unsigned char sign) {
  unsigned int i;

  for(i = 0; i < N; i++)
    a->coeffs[i] = CFLIP(a->coeffs[i], sign);

}

void poly_small(polysmall_t *a, const poly_t *b) {
  unsigned int i;
  int32_t bi;


  for(i = 0; i < N; i++) {
      bi = b->coeffs[i]; 
      a->coeffs[i] = CMUX(bi, bi-Q, 2*bi<Q);
  }

}

void poly_unsmall(poly_t *a, const polysmall_t *b) {
  unsigned int i;
  int32_t bi;


  for(i = 0; i < N; i++) {
      bi = b->coeffs[i]; 
      a->coeffs[i] = CMUX(bi, Q+bi, bi>=0);
  }

}

void poly_ntt(poly_t *a) {

  ntt(a->coeffs);

}

void poly_invntt_montgomery(poly_t *a) {

  invntt_frominvmont(a->coeffs);

}

void poly_pointwise_invmontgomery(poly_t *c, const poly_t *a, const poly_t *b) {
  unsigned int i;

  for(i = 0; i < N; ++i)
    c->coeffs[i] = montgomery_reduce((uint64_t)a->coeffs[i] * b->coeffs[i]);

}

void poly_pack4(unsigned char *r, const polysmall_t *restrict a,
	int32_t offset)
{
    uint32_t r0, r1;
    for(int i=0; i<N/2; i++) {
	r0 = a->coeffs[2*i+0] + offset;
	r1 = a->coeffs[2*i+1] + offset;

	r[i] = r0 | (r1 << 4);
    }	
}

int poly_unpack4(polysmall_t *restrict a, const unsigned char* r,
	int32_t offset)
{
    for(int i=0; i<N/2; i++) {
	a->coeffs[2*i+0] = (int32_t)(r[i] & 0xF) - offset;
	a->coeffs[2*i+1] = (int32_t)(r[i] >> 4)  - offset;
    }	
    return 0;
}

void poly_pack8(unsigned char *r, const polysmall_t *restrict a,
	int32_t offset)
{
    for(int i=0; i<N; i++) {
	r[i] = a->coeffs[i] + offset;
    }	
}

int poly_unpack8(polysmall_t *restrict a, const unsigned char* r,
	int32_t offset)
{
    for(int i=0; i<N; i++) {
	a->coeffs[i] = (int32_t)r[i] - offset;
    }	
    return 0;
}

void poly_pack16(unsigned char *r, const poly_t *restrict a)
{
    for(int i=0; i<N; i++) {
	r[2*i+0] = a->coeffs[i] & 0xFF;
	r[2*i+1] = a->coeffs[i] >> 8;
    }	
}

int poly_unpack16(poly_t *restrict a, const unsigned char* r)
{
    for(int i=0; i<N; i++) {
	a->coeffs[i] = r[2*i] | ((uint32_t)r[2*i+1] << 8);
    }	
    return 0;
}

void poly_rot(polysmall_t *restrict a, const polysmall_t *restrict b,
	int i)
{
    int j;
    for(int k=0; k<N; k++) {
	j = N-i+k;
	a->coeffs[k] = (j < N) ? -b->coeffs[j] : b->coeffs[j-N];
    }
}
