#include "packing.h"
#include "poly.h"
#include <string.h>

int pack_pk(unsigned char pk[CRYPTO_PUBLICKEYBYTES],
	const poly_t *restrict a_fft)
{
    poly_pack16(pk, a_fft);
    return 0;
}

int unpack_pk(poly_t *restrict a_fft,
	const unsigned char pk[CRYPTO_PUBLICKEYBYTES])
{
    poly_unpack16(a_fft, pk);
    return 0;
}

int pack_sk(unsigned char sk[CRYPTO_SECRETKEYBYTES],
	const unsigned char pk[CRYPTO_PUBLICKEYBYTES],
	const polysmall_t *restrict s1, const polysmall_t *restrict s2)
{
    memmove(sk, pk, CRYPTO_PUBLICKEYBYTES);
    poly_pack4(sk + CRYPTO_PUBLICKEYBYTES, s1, 4);
    poly_pack4(sk + CRYPTO_PUBLICKEYBYTES + PACK4_BYTES, s2, 4);
    return 0;
}

int unpack_sk(polysmall_t *restrict s1, polysmall_t *restrict s2,
	poly_t *restrict a_fft,
	const unsigned char sk[CRYPTO_SECRETKEYBYTES])
{
    poly_unpack16(a_fft, sk);
    poly_unpack4(s1, sk + CRYPTO_PUBLICKEYBYTES, 4);
    poly_unpack4(s2, sk + CRYPTO_PUBLICKEYBYTES + PACK4_BYTES, 4);
    return 0;
}

int pack_sk_fft(unsigned char sk[CRYPTO_SECRETKEYBYTES],
	const unsigned char pk[CRYPTO_PUBLICKEYBYTES],
	const poly_t *restrict s1_fft, const poly_t *restrict s2_fft)
{
    memmove(sk, pk, CRYPTO_PUBLICKEYBYTES);
    poly_pack16(sk + CRYPTO_PUBLICKEYBYTES, s1_fft);
    poly_pack16(sk + CRYPTO_PUBLICKEYBYTES + PACK16_BYTES, s2_fft);
    return 0;
}

int unpack_sk_fft(poly_t *restrict s1_fft, poly_t *restrict s2_fft,
	poly_t *restrict a_fft,
	const unsigned char sk[CRYPTO_SECRETKEYBYTES])
{
    poly_unpack16(a_fft, sk);
    poly_unpack16(s1_fft, sk + CRYPTO_PUBLICKEYBYTES);
    poly_unpack16(s2_fft, sk + CRYPTO_PUBLICKEYBYTES + PACK16_BYTES);

    /* could do consistency check on key sizes */
    return 0;
}

void pack_hpol(unsigned char *dest, const hpol_t *restrict c)
{
    for(int i=0; i<kappa; i++) {
	dest[2*i]   = c->indices[i] & 0xFF;
	dest[2*i+1] = c->indices[i] >> 8;
    }
}

int unpack_hpol(hpol_t *restrict c, const unsigned char *src)
{
    uint16_t ci;
    for(int i=0; i<kappa; i++) {
	ci = src[2*i] + (src[2*i+1] << 8);
	if(ci >= N)
	    return -1;
	c->indices[i] = ci;
    }
    return 0;
}

int pack_sig(unsigned char sm[CRYPTO_BYTES],
	const polysmall_t *restrict z1, const polysmall_t *restrict z2dag,
	const hpol_t *restrict c)
{
    poly_t zbig;
    poly_unsmall(&zbig, z1);
    poly_pack16(sm, &zbig);
    poly_pack8(sm + PACK16_BYTES, z2dag, modp);
    pack_hpol(sm + PACK16_BYTES + PACK8_BYTES, c);

    return 0;
}

int unpack_sig(polysmall_t *restrict z1, polysmall_t *restrict z2dag,
	hpol_t *restrict csig, const unsigned char *sm)
{
    int err = 0;

    poly_t zbig;
    poly_unpack16(&zbig, sm);
    poly_small(z1, &zbig);

    poly_unpack8(z2dag, sm + PACK16_BYTES, modp);
    err = unpack_hpol(csig, sm + PACK16_BYTES + PACK8_BYTES);

    return err;
}

int pack_sig_nocomp(unsigned char sm[CRYPTO_BYTES],
	const polysmall_t *restrict z1, const polysmall_t *restrict z2,
	const hpol_t *restrict c)
{
    poly_t zbig;
    poly_unsmall(&zbig, z1);
    poly_pack16(sm, &zbig);

    poly_unsmall(&zbig, z2);
    poly_pack16(sm + PACK16_BYTES, &zbig);

    pack_hpol(sm + 2*PACK16_BYTES, c);

    return 0;
}

int unpack_sig_nocomp(polysmall_t *restrict z1, polysmall_t *restrict z2,
	hpol_t *restrict csig, const unsigned char *sm)
{
    int err = 0;

    poly_t zbig;
    poly_unpack16(&zbig, sm);
    poly_small(z1, &zbig);

    poly_unpack16(&zbig, sm + PACK16_BYTES);
    poly_small(z2, &zbig);

    err = unpack_hpol(csig, sm + 2*PACK16_BYTES);

    return err;
}


