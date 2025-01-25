#include <arm_math.h>

#include <assert.h>

/**
 * @brief Conversion of quaternion to equivalent rotation matrix.
 * @param[in]       pInputQuaternions points to an array of normalized quaternions
 * @param[out]      pOutputRotations points to an array of 3x3 rotations (in row order)
 * @param[in]       nbQuaternions in the array
 * @return none.
 *
 * <b>Format of rotation matrix</b>
 * \par
 * The quaternion a + ib + jc + kd is converted into rotation matrix:
 *   a^2 + b^2 - c^2 - d^2                 2bc - 2ad                 2bd + 2ac
 *               2bc + 2ad     a^2 - b^2 + c^2 - d^2                 2cd - 2ab
 *               2bd - 2ac                 2cd + 2ab     a^2 - b^2 - c^2 + d^2
 *
 * Rotation matrix is saved in row order : R00 R01 R02 R10 R11 R12 R20 R21 R22
 */
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