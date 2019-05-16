#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

#define MONT 10952U // 2^32 % Q
#define QINV 4143984639U // -q^(-1) mod 2^32

uint32_t montgomery_reduce(uint64_t a);

uint32_t reduce32(uint32_t a);
uint32_t reduce2q(uint32_t a);
uint32_t reducep(uint32_t a);

uint32_t csubq(uint32_t a);
uint32_t csub2q(uint32_t a);
uint32_t csubp(uint32_t a);

uint32_t freeze(uint32_t a);

uint32_t inverse_modq(uint32_t a);

#endif
