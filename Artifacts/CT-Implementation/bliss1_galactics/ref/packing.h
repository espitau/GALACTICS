#ifndef PACKING_H
#define PACKING_H

#include <stdint.h>
#include "bliss.h"
#include "api.h"

#define PACK4_BYTES  256
#define PACK8_BYTES  512
#define PACK16_BYTES 1024

int unpack_sk(polysmall_t *restrict s1, polysmall_t *restrict s2,
	poly_t *restrict a_fft,
	const unsigned char sk[CRYPTO_SECRETKEYBYTES]);

int unpack_sk_fft(poly_t *restrict s1_fft, poly_t *restrict s2_fft,
	poly_t *restrict a_fft,
	const unsigned char sk[CRYPTO_SECRETKEYBYTES]);

int unpack_pk(poly_t *restrict a_fft,
	const unsigned char pk[CRYPTO_PUBLICKEYBYTES]);

int unpack_sig(polysmall_t *restrict z1, polysmall_t *restrict z2dag,
	hpol_t *restrict csig, const unsigned char *sm);

int pack_sig(unsigned char sm[CRYPTO_BYTES],
	const polysmall_t *restrict z1, const polysmall_t *restrict z2dag,
	const hpol_t *restrict csig);

int pack_pk(unsigned char pk[CRYPTO_PUBLICKEYBYTES],
	const poly_t *restrict a_fft);

int pack_sk(unsigned char sk[CRYPTO_SECRETKEYBYTES],
	const unsigned char pk[CRYPTO_PUBLICKEYBYTES],
	const polysmall_t *restrict s1, const polysmall_t *restrict s2);

int pack_sk_fft(unsigned char sk[CRYPTO_SECRETKEYBYTES],
	const unsigned char pk[CRYPTO_PUBLICKEYBYTES],
	const poly_t *restrict s1_fft, const poly_t *restrict s2_fft);

int pack_sig_nocomp(unsigned char sm[CRYPTO_BYTES],
	const polysmall_t *restrict z1, const polysmall_t *restrict z2,
	const hpol_t *restrict c);

int unpack_sig_nocomp(polysmall_t *restrict z1, polysmall_t *restrict z2,
	hpol_t *restrict csig, const unsigned char *sm);
#endif
