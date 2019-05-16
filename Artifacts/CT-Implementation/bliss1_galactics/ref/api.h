#define COMPRESS		1
#define MULTC_FFT		1
#define CRYPTO_PUBLICKEYBYTES 	1024  // == 2*N

#if MULTC_FFT
#define CRYPTO_SECRETKEYBYTES 	3072  // == 6*N
#else
#define CRYPTO_SECRETKEYBYTES 	1536  // == 2*N + N
#endif

#if COMPRESS
#define CRYPTO_BYTES		1582  // == 2*N + N + 2*kappa
#else
#define CRYPTO_BYTES		2094  // == 2*N + 2*N + 2*kappa
#endif

#define CRYPTO_DETERMINISTIC	0

