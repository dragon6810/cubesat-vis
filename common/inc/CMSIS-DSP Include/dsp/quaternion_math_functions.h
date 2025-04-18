#ifndef _quaternion_math_functions_h
#define _quaternion_math_functions_h

#include <CMSIS-DSP Include/arm_math_memory.h>
#include <CMSIS-DSP Include/arm_math_types.h>

/**
  @brief         Floating-point quaternion Norm.
  @param[in]     pInputQuaternions       points to the input vector of quaternions
  @param[out]    pNorms                  points to the output vector of norms
  @param[in]     nbQuaternions           number of quaternions in each vector
  @return        none
 */
void arm_quaternion_norm_f32(const float32_t *pInputQuaternions, 
    float32_t *pNorms,
    uint32_t nbQuaternions);

/**
  @brief         Floating-point quaternion inverse.
  @param[in]     pInputQuaternions            points to the input vector of quaternions
  @param[out]    pInverseQuaternions          points to the output vector of inverse quaternions
  @param[in]     nbQuaternions                number of quaternions in each vector
  @return        none
 */
void arm_quaternion_inverse_f32(const float32_t *pInputQuaternions, 
    float32_t *pInverseQuaternions, 
    uint32_t nbQuaternions);

/**
  @brief         Floating-point quaternion conjugates.
  @param[in]     pInputQuaternions            points to the input vector of quaternions
  @param[out]    pConjugateQuaternions        points to the output vector of conjugate quaternions
  @param[in]     nbQuaternions                number of quaternions in each vector
  @return        none
 */
void arm_quaternion_conjugate_f32(const float32_t *inputQuaternions, 
    float32_t *pConjugateQuaternions, 
    uint32_t nbQuaternions);

/**
  @brief         Floating-point normalization of quaternions.
  @param[in]     pInputQuaternions            points to the input vector of quaternions
  @param[out]    pNormalizedQuaternions       points to the output vector of normalized quaternions
  @param[in]     nbQuaternions                number of quaternions in each vector
  @return        none
 */
void arm_quaternion_normalize_f32(const float32_t *inputQuaternions, 
    float32_t *pNormalizedQuaternions, 
    uint32_t nbQuaternions);

/**
  @brief         Floating-point product of two quaternions.
  @param[in]     qa       First quaternion
  @param[in]     qb       Second quaternion
  @param[out]    r        Product of two quaternions
  @return        none
 */
void arm_quaternion_product_single_f32(const float32_t *qa, 
    const float32_t *qb, 
    float32_t *r);

/**
  @brief         Floating-point elementwise product two quaternions.
  @param[in]     qa                  First array of quaternions
  @param[in]     qb                  Second array of quaternions
  @param[out]    r                   Elementwise product of quaternions
  @param[in]     nbQuaternions       Number of quaternions in the array
  @return        none
 */
void arm_quaternion_product_f32(const float32_t *qa, 
    const float32_t *qb, 
    float32_t *r,
    uint32_t nbQuaternions);

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
    uint32_t nbQuaternions);

/**
 * @brief Conversion of a rotation matrix to equivalent quaternion.
 * @param[in]       pInputRotations points to an array 3x3 rotation matrix (in row order)
 * @param[out]      pOutputQuaternions points to an array of quaternions
 * @param[in]       nbQuaternions in the array
 * @return none.
*/
void arm_rotation2quaternion_f32(const float32_t *pInputRotations, 
    float32_t *pOutputQuaternions,  
    uint32_t nbQuaternions);

#endif