#include <VecUtils.h>

#include <stdio.h>
#include <assert.h>
#include <math.h>

void Vec_RotateSpher(Vec3D_t* V, float theta, float phi, Vec3D_t* result)
{
    Vec3D_t res;
    Vec3D_t temp;
    float initialangle;
	
    assert(V);
	
    initialangle = -atan2f(V->Y, V->X);
    
    #if 1
    // Rotate res onto the XZ plane
    temp = *V;
    res.X = temp.X * cosf(initialangle) - temp.Y * sinf(initialangle);
    res.Y = temp.X * sinf(initialangle) + temp.Y * cosf(initialangle);
    res.Z = temp.Z;

    // Rotate res along the XZ plane by theta (up/down)
    temp = res;
    res.X = temp.X * cosf(theta) - temp.Z * sinf(theta);
    res.Z = temp.X * sinf(theta) + temp.Z * cosf(theta);

    // Rotate it back and by phi
    initialangle = -initialangle + phi;
    temp = res;
    res.X = temp.X * cosf(initialangle) - temp.Y * sinf(initialangle);
    res.Y = temp.X * sinf(initialangle) + temp.Y * cosf(initialangle);
    #else
        temp = *V;
        res.X = temp.X * cosf(phi) - temp.Y * sinf(phi);
        res.Y = temp.X * sinf(phi) + temp.Y * cosf(phi);
        res.Z = temp.Z;
    #endif

    if(result)
        *result = res;
    else
        *V = res;
}
