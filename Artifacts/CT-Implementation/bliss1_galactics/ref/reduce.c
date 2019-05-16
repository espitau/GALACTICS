#include <inttypes.h>
#include "bliss.h"
#include "reduce.h"

#include <stdio.h>

uint32_t inverse_modq(uint32_t a) {
    /* can easily be converted to constant-time
     * (but we do not claim constant time key generation, so not needed)
     */
    uint32_t u = Q, v = a, r = 0, s = 1, k = 0;

    while( v > 1 ) {
	if( (u&1) == 0 ) {
	    u >>= 1;
	    s <<= 1;
	}
	else if( (v&1) == 0 ) {
	    v >>= 1;
	    r <<= 1;
	} else if( u > v ) {
	    u   = (u-v) >> 1;
	    r  += s;
	    s <<= 1;
	} else {
	    v   = (v-u) >> 1;
	    s  += r;
	    r <<= 1;
	}
	k++;
    }

    s = montgomery_reduce((uint64_t)s << (32-k));
    return s;
}
uint32_t montgomery_reduce(uint64_t a) {
  uint64_t t;

  t = a * QINV;
  t &= (1ULL << 32) - 1;
  t *= Q;
  t = a + t;
  t >>= 32;
  return t;
}

uint32_t reduce32(uint32_t a) {
  uint32_t t;
  t = ((uint64_t)2863078533LL * a) >> 45;
  return a - Q*t;
}

uint32_t reduce2q(uint32_t a) {
  uint32_t t;
  t = ((uint64_t)2863078533ULL * a) >> 46;
  return a - 2*Q*t;
}

uint32_t reducep(uint32_t a) {
  uint32_t t;
  t = ((uint64_t)2863311531ULL * a) >> 36;
  return a - modp*t;
}

uint32_t csubq(uint32_t a) {
  a -= Q;
  a += ((int32_t)a >> 31) & Q;
  return a;
}

uint32_t csub2q(uint32_t a) {
  a -= 2*Q;
  a += ((int32_t)a >> 31) & (2*Q);
  return a;
}

uint32_t csubp(uint32_t a) {
  a -= modp;
  a += ((int32_t)a >> 31) & modp;
  return a;
}

uint32_t freeze(uint32_t a) {
  a = reduce32(a);
  return a;
}
