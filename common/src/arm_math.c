#include <arm_math.h>

arm_status arm_sqrt_f32(const float in, float* out)
{
    *out = sqrtf(in);
    
    if(in < 0)
        return ARM_MATH_NANINF;

    return ARM_MATH_SUCCESS; 
}

float arm_sin_f32(const float x)
{
    return sinf(x);
}

float arm_cos_f32(const float x)
{
    return cosf(x);
}

arm_status arm_atan2_f32(float y, float x, float* result)
{
    *result = atan2f(y, x);
    return ARM_MATH_SUCCESS;
}
