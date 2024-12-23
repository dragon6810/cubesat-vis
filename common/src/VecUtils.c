#include <VecUtils.h>

#include <assert.h>

void Vec_RotateSpher(Vec3D_t* V, float theta, float phi, Vec3D_t* result)
{
    Vec3D_t res;
	
    assert(V);
	
	float sin_theta = arm_sin_f32(theta); 
    float cos_theta = arm_cos_f32(theta);

	float sin_phi = arm_sin_f32(phi);
    float cos_phi = arm_cos_f32(phi);
	
    res.X = V->X * cos_phi * cos_theta
		+ V->Z * cos_phi * sin_theta
		- V->Y * sin_phi;
	
    res.Y = V->X * sin_phi * cos_theta
		+ V->Z * sin_phi * sin_theta
		+ V->Y * cos_phi;
	
    res.Z = V->Z * cos_theta
		- V->X * sin_theta;

    memcpy(result != NULL ? result : V, &res, sizeof(Vec3D_t));
}
