#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "bliss.h"
#include "fips202.h"

void poly_add(polysmall_t *restrict a, const polysmall_t *restrict b);

void poly_ntt(poly_t *a);
void poly_invntt_montgomery(poly_t *a);
void poly_pointwise_invmontgomery(poly_t *c, const poly_t *a, const poly_t *b);

void poly_csubq(poly_t *a);
void poly_reduce(poly_t *a);
void poly_freeze(poly_t *a);

void poly_sample_sparse(polysmall_t *a);

void poly_small(polysmall_t *a, const poly_t *b);
void poly_unsmall(poly_t *a, const polysmall_t *b);
void poly_sign_flip(polysmall_t *a, unsigned char bit);

void poly_pack4(unsigned char* dest, const polysmall_t *restrict a,
	int32_t offset);
void poly_pack8(unsigned char* dest, const polysmall_t *restrict a,
	int32_t offset);
void poly_pack16(unsigned char* dest, const poly_t *restrict a);

int poly_unpack4(polysmall_t *restrict a, const unsigned char* src,
	int32_t offset);
int poly_unpack8(polysmall_t *restrict a, const unsigned char* src,
	int32_t offset);
int poly_unpack16(poly_t *restrict a, const unsigned char* src);

void poly_rot(polysmall_t *restrict a, const polysmall_t *restrict b, int i);
#endif
