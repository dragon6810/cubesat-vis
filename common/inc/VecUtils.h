#ifndef _vecutils_h
#define _vecutils_h

#include <User_types.h>

#ifdef __cplusplus
extern "C"
{
#endif

void Vec_RotateSpher(Vec3D_t* V, float theta, float phi, Vec3D_t* result);

#ifdef __cplusplus
}
#endif

#endif
