#include <VecUtils.h>

#include <assert.h>
#include <math.h>

void Vec_RotateSpher(Vec3D_t* V, float theta, float phi, Vec3D_t* result)
{
    Vec3D_t res;
    Vec3D_t before;
    float initialangle;
	
    assert(V);
	
    initialangle = -atan2f(V->Y, V->X);
    
    // Rotate res onto the XZ plane
    res.X = V->X * cosf(initialangle) - V->Y * sinf(initialangle);
    res.Y = V->X * sinf(initialangle) + V->Y * cosf(initialangle);

    // Rotate res along the XZ plane by theta (up/down)
    before = res;
    res.X = before.X * cosf(theta) - before.Z * sinf(theta);
    res.Z = before.X * sinf(theta) + before.Z * cosf(theta);

    // Rotate it back and by phi
    initialangle = -initialangle + phi;
    before = res;
    res.X = before.X * cosf(initialangle) - before.Y * sinf(initialangle);
    res.Y = before.X * sinf(initialangle) + before.Y * cosf(initialangle);

    if(result)
        *result = res;
    else
        *V = res;
}
