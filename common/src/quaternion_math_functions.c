#include <arm_math.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void arm_quaternion2rotation_f32(const float32_t *pInputQuaternions, 
    float32_t *pOutputRotations, 
    uint32_t nbQuaternions)
{
    int i;

    float32_t *q, *mat;
    float32_t a, b, c, d;

    if(!nbQuaternions)
        return;

    assert(pInputQuaternions);
    assert(pOutputRotations);

    for(i=0, q=pInputQuaternions, mat=pOutputRotations; i<nbQuaternions; i++, q+=4, mat+=9)
    {
        a = q[0]; b = q[1]; c = q[2]; d = q[3];

        // Row 1
        mat[0] = a*a + b*b - c*c - d*d;
        mat[1] = 2.0*b*c - 2.0*a*d;
        mat[2] = 2.0*b*d + 2.0*a*c;

        // Row 2
        mat[3] = 2.0*b*c + 2.0*a*d;
        mat[4] = a*a - b*b + c*c - d*d;
        mat[5] = 2.0*c*d - 2.0*a*b;

        // Row 3
        mat[6] = 2.0*b*d - 2.0*a*c;
        mat[7] = 2.0*c*d + 2.0*a*b;
        mat[8] = a*a - b*b - c*c + d*d;
    }
}

void arm_rotation2quaternion_f32(const float32_t *pInputRotations, 
    float32_t *pOutputQuaternions,  
    uint32_t nbQuaternions)
{
    printf("arm_rotation2quaternion_f32 not implemented yet.\n");
    abort();
}