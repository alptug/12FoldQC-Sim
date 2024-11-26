#pragma once
/**
 * @brief double precision SIMD oriented Fast Mersenne Twister(dSFMT)
 * pseudorandom number generator based on IEEE 754 format.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and
 * Hiroshima University. All rights reserved.
 * Copyright (C) 2012 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The new BSD License is applied to this software.
 * see LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t
 #define PRIu64 "llu"
 #define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#ifndef DSFMT_H
#define DSFMT_H
#if defined(__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#define DSFMT_MEXP 19937
/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/* Mersenne Exponent. The period of the sequence
 *  is a multiple of 2^DSFMT_MEXP-1.
 * #define DSFMT_MEXP 19937 */
/** DSFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define DSFMT_N ((DSFMT_MEXP - 128) / 104 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define DSFMT_N32 (DSFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

#if !defined(DSFMT_BIG_ENDIAN)
#  if defined(__BYTE_ORDER) && defined(__BIG_ENDIAN)
#    if __BYTE_ORDER == __BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(_BYTE_ORDER) && defined(_BIG_ENDIAN)
#    if _BYTE_ORDER == _BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BYTE_ORDER__) && defined(__BIG_ENDIAN__)
#    if __BYTE_ORDER__ == __BIG_ENDIAN__
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(BYTE_ORDER) && defined(BIG_ENDIAN)
#    if BYTE_ORDER == BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BIG_ENDIAN) || defined(_BIG_ENDIAN) \
    || defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN)
#      define DSFMT_BIG_ENDIAN 1
#  endif
#endif

#if defined(DSFMT_BIG_ENDIAN) && defined(__amd64)
#  undef DSFMT_BIG_ENDIAN
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#  include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
#  if !defined(DSFMT_UINT32_DEFINED) && !defined(SFMT_UINT32_DEFINED)
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#    ifndef UINT64_C
#      define UINT64_C(v) (v ## ui64)
#    endif
#    define DSFMT_UINT32_DEFINED
#    if !defined(inline) && !defined(__cplusplus)
#      define inline __inline
#    endif
#  endif
#else
#  include <inttypes.h>
#  if !defined(inline) && !defined(__cplusplus)
#    if defined(__GNUC__)
#      define inline __inline__
#    else
#      define inline
#    endif
#  endif
#endif

#ifndef PRIu64
#  if defined(_MSC_VER) || defined(__BORLANDC__)
#    define PRIu64 "I64u"
#    define PRIx64 "I64x"
#  else
#    define PRIu64 "llu"
#    define PRIx64 "llx"
#  endif
#endif

#ifndef UINT64_C
#  define UINT64_C(v) (v ## ULL)
#endif

/*----------------------
  the parameters of DSFMT
  following definitions are in dSFMT-paramsXXXX.h file.
  ----------------------*/
/** the pick up position of the array.
#define DSFMT_POS1 122
*/

/** the parameter of shift left as four 32-bit registers.
#define DSFMT_SL1 18
 */

/** the parameter of shift right as four 32-bit registers.
#define DSFMT_SR1 12
*/

/** A bitmask, used in the recursion.  These parameters are introduced
 * to break symmetry of SIMD.
#define DSFMT_MSK1 (uint64_t)0xdfffffefULL
#define DSFMT_MSK2 (uint64_t)0xddfecb7fULL
*/

/** These definitions are part of a 128-bit period certification vector.
#define DSFMT_PCV1    UINT64_C(0x00000001)
#define DSFMT_PCV2    UINT64_C(0x00000000)
*/

#define DSFMT_LOW_MASK  UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_SR    12

/* for sse2 */
#if defined(HAVE_SSE2)
  #define SSE2_SHUFF 0x1b
#elif defined(HAVE_ALTIVEC)
  #if defined(__APPLE__)  /* For OSX */
    #define ALTI_SR (vector unsigned char)(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
    #define ALTI_SR_PERM \
        (vector unsigned char)(15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14)
    #define ALTI_SR_MSK \
        (vector unsigned int)(0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU)
    #define ALTI_PERM \
        (vector unsigned char)(12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3)
  #else
    #define ALTI_SR      {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}
    #define ALTI_SR_PERM {15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14}
    #define ALTI_SR_MSK  {0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU}
    #define ALTI_PERM    {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3}
  #endif
#endif
/* #define DSFMT_N    191 */
/* #define DSFMT_MAXDEGREE    19992 */
#define DSFMT_POS1    117
#define DSFMT_SL1    19
#define DSFMT_MSK1    UINT64_C(0x000ffafffffffb3f)
#define DSFMT_MSK2    UINT64_C(0x000ffdfffc90fffd)
#define DSFMT_MSK32_1    0x000ffaffU
#define DSFMT_MSK32_2    0xfffffb3fU
#define DSFMT_MSK32_3    0x000ffdffU
#define DSFMT_MSK32_4    0xfc90fffdU
#define DSFMT_FIX1    UINT64_C(0x90014964b32f4329)
#define DSFMT_FIX2    UINT64_C(0x3b8d12ac548a7c7a)
#define DSFMT_PCV1    UINT64_C(0x3d84e1ac0dc82880)
#define DSFMT_PCV2    UINT64_C(0x0000000000000001)
#define DSFMT_IDSTR    "dSFMT2-19937:117-19:ffafffffffb3f-ffdfffc90fffd"


/* PARAMETERS FOR ALTIVEC */
#if defined(__APPLE__)    /* For OSX */
    #define ALTI_SL1     (vector unsigned char)(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
    #define ALTI_SL1_PERM \
    (vector unsigned char)(2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1)
    #define ALTI_SL1_MSK \
    (vector unsigned int)(0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U)
    #define ALTI_MSK    (vector unsigned int)(DSFMT_MSK32_1, \
            DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4)
#else    /* For OTHER OSs(Linux?) */
    #define ALTI_SL1     {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}
    #define ALTI_SL1_PERM \
    {2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1}
    #define ALTI_SL1_MSK \
    {0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U}
    #define ALTI_MSK \
    {DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4}
#endif


/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
#  if !defined(__APPLE__)
#    include <altivec.h>
#  endif
/** 128-bit data structure */
union W128_T {
    vector unsigned int s;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};

#elif defined(HAVE_SSE2)
#  include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
    __m128i si;
    __m128d sd;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#else  /* standard C */
/** 128-bit data structure */
union W128_T {
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/** the 128-bit internal state array */
struct DSFMT_T {
    w128_t status[DSFMT_N + 1];
    int idx;
};
typedef struct DSFMT_T dsfmt_t;
/**
 * @file dSFMT-common.h
 */
#if defined(HAVE_SSE2)
#  include <emmintrin.h>
union X128I_T {
    uint64_t u[2];
    __m128i  i128;
};
union X128D_T {
    double d[2];
    __m128d d128;
};
/** mask data for sse2 */
static const union X128I_T sse2_param_mask = {{DSFMT_MSK1, DSFMT_MSK2}};
#endif

#if defined(HAVE_ALTIVEC)
inline static void do_recursion(w128_t *r, w128_t *a, w128_t * b,
                w128_t *lung) {
    const vector unsigned char sl1 = ALTI_SL1;
    const vector unsigned char sl1_perm = ALTI_SL1_PERM;
    const vector unsigned int sl1_msk = ALTI_SL1_MSK;
    const vector unsigned char sr1 = ALTI_SR;
    const vector unsigned char sr1_perm = ALTI_SR_PERM;
    const vector unsigned int sr1_msk = ALTI_SR_MSK;
    const vector unsigned char perm = ALTI_PERM;
    const vector unsigned int msk1 = ALTI_MSK;
    vector unsigned int w, x, y, z;

    z = a->s;
    w = lung->s;
    x = vec_perm(w, (vector unsigned int)perm, perm);
    y = vec_perm(z, (vector unsigned int)sl1_perm, sl1_perm);
    y = vec_sll(y, sl1);
    y = vec_and(y, sl1_msk);
    w = vec_xor(x, b->s);
    w = vec_xor(w, y);
    x = vec_perm(w, (vector unsigned int)sr1_perm, sr1_perm);
    x = vec_srl(x, sr1);
    x = vec_and(x, sr1_msk);
    y = vec_and(w, msk1);
    z = vec_xor(z, y);
    r->s = vec_xor(z, x);
    lung->s = w;
}
#elif defined(HAVE_SSE2)
/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param d a 128-bit part of the internal state array (I/O)
 */
inline static void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *u) {
    __m128i v, w, x, y, z;

    x = a->si;
    z = _mm_slli_epi64(x, DSFMT_SL1);
    y = _mm_shuffle_epi32(u->si, SSE2_SHUFF);
    z = _mm_xor_si128(z, b->si);
    y = _mm_xor_si128(y, z);

    v = _mm_srli_epi64(y, DSFMT_SR);
    w = _mm_and_si128(y, sse2_param_mask.i128);
    v = _mm_xor_si128(v, x);
    v = _mm_xor_si128(v, w);
    r->si = v;
    u->si = y;
}
#else
/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param lung a 128-bit part of the internal state array (I/O)
 */
inline static void do_recursion(w128_t *r, w128_t *a, w128_t * b,
                w128_t *lung) {
    uint64_t t0, t1, L0, L1;

    t0 = a->u[0];
    t1 = a->u[1];
    L0 = lung->u[0];
    L1 = lung->u[1];
    lung->u[0] = (t0 << DSFMT_SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b->u[0];
    lung->u[1] = (t1 << DSFMT_SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b->u[1];
    r->u[0] = (lung->u[0] >> DSFMT_SR) ^ (lung->u[0] & DSFMT_MSK1) ^ t0;
    r->u[1] = (lung->u[1] >> DSFMT_SR) ^ (lung->u[1] & DSFMT_MSK2) ^ t1;
}
#endif

/** dsfmt mexp for check */
extern const int dsfmt_global_mexp;

void dsfmt_gen_rand_all(dsfmt_t *dsfmt);
void dsfmt_chk_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed, int mexp);

#if defined(__GNUC__)
#  define DSFMT_PRE_INLINE inline static
#  define DSFMT_PST_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#  define DSFMT_PRE_INLINE __forceinline static
#  define DSFMT_PST_INLINE
#else
#  define DSFMT_PRE_INLINE inline static
#  define DSFMT_PST_INLINE
#endif
DSFMT_PRE_INLINE uint32_t dsfmt_genrand_uint32(dsfmt_t *dsfmt) DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_close_open(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_close(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE double dsfmt_genrand_open_open(dsfmt_t *dsfmt)
    DSFMT_PST_INLINE;
DSFMT_PRE_INLINE void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed)
    DSFMT_PST_INLINE;

/**
 * This function generates and returns unsigned 32-bit integer.
 * This is slower than SFMT, only for convenience usage.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static uint32_t dsfmt_genrand_uint32(dsfmt_t *dsfmt) {
    uint32_t r;
    uint64_t *psfmt64 = &dsfmt->status[0].u[0];

    if (dsfmt->idx >= DSFMT_N64) {
        dsfmt_gen_rand_all(dsfmt);
        dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++] & 0xffffffffU;
    return r;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [1, 2).  This is
 * the primitive and faster than generating numbers in other ranges.
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close1_open2(dsfmt_t *dsfmt) {
    double r;
    double *psfmt64 = &dsfmt->status[0].d[0];

    if (dsfmt->idx >= DSFMT_N64) {
        dsfmt_gen_rand_all(dsfmt);
        dsfmt->idx = 0;
    }
    r = psfmt64[dsfmt->idx++];
    return r;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range [0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_close_open(dsfmt_t *dsfmt) {
    return dsfmt_genrand_close1_open2(dsfmt) - 1.0;
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1].
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_close(dsfmt_t *dsfmt) {
    return 2.0 - dsfmt_genrand_close1_open2(dsfmt);
}

/**
 * This function generates and returns double precision pseudorandom
 * number which distributes uniformly in the range (0, 1).
 * dsfmt_init_gen_rand() or dsfmt_init_by_array() must be called
 * before this function.
 * @param dsfmt dsfmt internal state date
 * @return double precision floating point pseudorandom number
 */
inline static double dsfmt_genrand_open_open(dsfmt_t *dsfmt) {
    double *dsfmt64 = &dsfmt->status[0].d[0];
    union {
        double d;
        uint64_t u;
    } r;

    if (dsfmt->idx >= DSFMT_N64) {
        dsfmt_gen_rand_all(dsfmt);
        dsfmt->idx = 0;
    }
    r.d = dsfmt64[dsfmt->idx++];
    r.u |= 1;
    return r.d - 1.0;
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 */
inline static void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed) {
    dsfmt_chk_init_gen_rand(dsfmt, seed, DSFMT_MEXP);
}

/* @file dSFMT.c */
/** dsfmt mexp for check */
static const int dsfmt_mexp = DSFMT_MEXP;

/*----------------
  STATIC FUNCTIONS
  ----------------*/
inline static uint32_t ini_func1(uint32_t x);
inline static uint32_t ini_func2(uint32_t x);
inline static void gen_rand_array_c1o2(dsfmt_t *dsfmt, w128_t *array,
                       int size);
inline static void gen_rand_array_c0o1(dsfmt_t *dsfmt, w128_t *array,
                       int size);
inline static void gen_rand_array_o0c1(dsfmt_t *dsfmt, w128_t *array,
                       int size);
inline static void gen_rand_array_o0o1(dsfmt_t *dsfmt, w128_t *array,
                       int size);
inline static int idxof(int i);
static void initial_mask(dsfmt_t *dsfmt);
static void period_certification(dsfmt_t *dsfmt);

#if defined(HAVE_SSE2)
/** 1 in 64bit for sse2 */
static const union X128I_T sse2_int_one = {{1, 1}};
/** 2.0 double for sse2 */
static const union X128D_T sse2_double_two = {{2.0, 2.0}};
/** -1.0 double for sse2 */
static const union X128D_T sse2_double_m_one = {{-1.0, -1.0}};
#endif

/**
 * This function simulate a 32-bit array index overlapped to 64-bit
 * array of LITTLE ENDIAN in BIG ENDIAN machine.
 */
#if defined(DSFMT_BIG_ENDIAN)
inline static int idxof(int i) {
    return i ^ 1;
}
#else
inline static int idxof(int i) {
    return i;
}
#endif

#if defined(HAVE_SSE2)
/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range [0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_c0o1(w128_t *w) {
    w->sd = _mm_add_pd(w->sd, sse2_double_m_one.d128);
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1].
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0c1(w128_t *w) {
    w->sd = _mm_sub_pd(sse2_double_two.d128, w->sd);
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0o1(w128_t *w) {
    w->si = _mm_or_si128(w->si, sse2_int_one.i128);
    w->sd = _mm_add_pd(w->sd, sse2_double_m_one.d128);
}
#else /* standard C and altivec */
/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range [0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_c0o1(w128_t *w) {
    w->d[0] -= 1.0;
    w->d[1] -= 1.0;
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1].
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0c1(w128_t *w) {
    w->d[0] = 2.0 - w->d[0];
    w->d[1] = 2.0 - w->d[1];
}

/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range (0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_o0o1(w128_t *w) {
    w->u[0] |= 1;
    w->u[1] |= 1;
    w->d[0] -= 1.0;
    w->d[1] -= 1.0;
}
#endif

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_c1o2(dsfmt_t *dsfmt, w128_t *array,
                       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
         &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
    dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    dsfmt->status[j] = array[i];
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_c0o1(dsfmt_t *dsfmt, w128_t *array,
                       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
         &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    convert_c0o1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
    dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    dsfmt->status[j] = array[i];
    convert_c0o1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
    convert_c0o1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_o0o1(dsfmt_t *dsfmt, w128_t *array,
                       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
         &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    convert_o0o1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
    dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    dsfmt->status[j] = array[i];
    convert_o0o1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
    convert_o0o1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_o0c1(dsfmt_t *dsfmt, w128_t *array,
                       int size) {
    int i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1],
         &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
    do_recursion(&array[i], &dsfmt->status[i],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    convert_o0c1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
    dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
    do_recursion(&array[i], &array[i - DSFMT_N],
             &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    dsfmt->status[j] = array[i];
    convert_o0c1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
    convert_o0c1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function represents a function used in the initialization
 * by init_by_array
 * @param x 32-bit integer
 * @return 32-bit integer
 */
static uint32_t ini_func1(uint32_t x) {
    return (x ^ (x >> 27)) * (uint32_t)1664525UL;
}

/**
 * This function represents a function used in the initialization
 * by init_by_array
 * @param x 32-bit integer
 * @return 32-bit integer
 */
static uint32_t ini_func2(uint32_t x) {
    return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
}

/**
 * This function initializes the internal state array to fit the IEEE
 * 754 format.
 * @param dsfmt dsfmt state vector.
 */
static void initial_mask(dsfmt_t *dsfmt) {
    int i;
    uint64_t *psfmt;

    psfmt = &dsfmt->status[0].u[0];
    for (i = 0; i < DSFMT_N * 2; i++) {
        psfmt[i] = (psfmt[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;
    }
}

/**
 * This function certificate the period of 2^{SFMT_MEXP}-1.
 * @param dsfmt dsfmt state vector.
 */
static void period_certification(dsfmt_t *dsfmt) {
    uint64_t pcv[2] = {DSFMT_PCV1, DSFMT_PCV2};
    uint64_t tmp[2];
    uint64_t inner;
    int i;
#if (DSFMT_PCV2 & 1) != 1
    int j;
    uint64_t work;
#endif

    tmp[0] = (dsfmt->status[DSFMT_N].u[0] ^ DSFMT_FIX1);
    tmp[1] = (dsfmt->status[DSFMT_N].u[1] ^ DSFMT_FIX2);

    inner = tmp[0] & pcv[0];
    inner ^= tmp[1] & pcv[1];
    for (i = 32; i > 0; i >>= 1) {
        inner ^= inner >> i;
    }
    inner &= 1;
    /* check OK */
    if (inner == 1) {
    return;
    }
    /* check NG, and modification */
#if (DSFMT_PCV2 & 1) == 1
    dsfmt->status[DSFMT_N].u[1] ^= 1;
#else
    for (i = 1; i >= 0; i--) {
    work = 1;
    for (j = 0; j < 64; j++) {
        if ((work & pcv[i]) != 0) {
        dsfmt->status[DSFMT_N].u[i] ^= work;
        return;
        }
        work = work << 1;
    }
    }
#endif
    return;
}

/*----------------
  PUBLIC FUNCTIONS
  ----------------*/


/**
 * This function fills the internal state array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 */
void dsfmt_gen_rand_all(dsfmt_t *dsfmt) {
    int i;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&dsfmt->status[0], &dsfmt->status[0],
         &dsfmt->status[DSFMT_POS1], &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
    do_recursion(&dsfmt->status[i], &dsfmt->status[i],
             &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
    do_recursion(&dsfmt->status[i], &dsfmt->status[i],
             &dsfmt->status[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    dsfmt->status[DSFMT_N] = lung;
}


#if defined(__INTEL_COMPILER)
#  pragma warning(disable:981)
#endif
/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 * @param mexp caller's mersenne expornent
 */
void dsfmt_chk_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed, int mexp) {
    int i;
    uint32_t *psfmt;

    /* make sure caller program is compiled with the same MEXP */
    if (mexp != dsfmt_mexp) {
    fprintf(stderr, "DSFMT_MEXP doesn't match with dSFMT.c\n");
    exit(1);
    }
    psfmt = &dsfmt->status[0].u32[0];
    psfmt[idxof(0)] = seed;
    for (i = 1; i < (DSFMT_N + 1) * 4; i++) {
        psfmt[idxof(i)] = 1812433253UL
        * (psfmt[idxof(i - 1)] ^ (psfmt[idxof(i - 1)] >> 30)) + i;
    }
    initial_mask(dsfmt);
    period_certification(dsfmt);
    dsfmt->idx = DSFMT_N64;
}

#if defined(__INTEL_COMPILER)
#  pragma warning(default:981)
#endif

/* Global PRNG */
dsfmt_t dsfmt_global_data;

static inline void dsfmt_seed(uint32_t seed){
    dsfmt_init_gen_rand(&dsfmt_global_data, seed);
}

static inline double dsfmt_genrand(void){
    return dsfmt_genrand_close_open(&dsfmt_global_data);
}

#if defined(__cplusplus)
}
#endif

#endif /* DSFMT_H */
