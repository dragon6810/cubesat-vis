#ifndef _arm_math_types_h
#define _arm_math_types_h

#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef float float32_t;
typedef double float64_t;

typedef enum
{
    ARM_MATH_SUCCESS =                0,
    ARM_MATH_ARGUMENT_ERROR =        -1,
    ARM_MATH_LENGTH_ERROR =          -2,
    ARM_MATH_SIZE_MISMATCH =         -3,
    ARM_MATH_NANINF =                -4,
    ARM_MATH_SINGULAR =              -5,
    ARM_MATH_TEST_FAILURE =          -6,
    ARM_MATH_DECOMPOSITION_FAILURE = -7,
} arm_status;

#ifdef __cplusplus
}
#endif

#endif
