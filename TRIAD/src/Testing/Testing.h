#ifndef _testing_h
#define _testing_h

#include <User_types.h>

void Testing_TestVectorImpl(Vec3D_t expect, Vec3D_t in, Vec3D_t out, const char* test, const char* file, int line);

#define Testing_TestVector(exp, in, out, test) (Testing_TestVectorImpl(exp, in, out, test, __FILE__, __LINE__))

#endif