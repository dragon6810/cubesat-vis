#ifndef _basic_math_functions_h
#define _basic_math_functions_h

#include <CMSIS-DSP Include/arm_math_types.h>

#ifdef   __cplusplus
extern "C"
{
#endif

inline arm_status arm_sqrt_f32(const float in, float* out) __attribute__((used));
inline float arm_sin_f32(const float x) __attribute__((used));
inline float arm_cos_f32(const float x) __attribute__((used));
inline arm_status arm_atan2_f32(float y, float x, float* result) __attribute__((used));

#ifdef   __cplusplus
}
#endif

#endif
