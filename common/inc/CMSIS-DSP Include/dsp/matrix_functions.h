#ifndef _matrix_functions_h
#define _matrix_functions_h

#include <CMSIS-DSP Include/arm_math_memory.h>
#include <CMSIS-DSP Include/arm_math_types.h>

#ifdef   __cplusplus
extern "C"
{
#endif

/**
 * @brief Instance structure for the floating-point matrix structure.
 */
typedef struct
{
    uint16_t numRows;     /**< number of rows of the matrix.     */
    uint16_t numCols;     /**< number of columns of the matrix.  */
    float32_t *pData;     /**< points to the data of the matrix. */
} arm_matrix_instance_f32;

/**
 * @brief Floating-point matrix addition.
 * @param[in]  pSrcA  points to the first input matrix structure
 * @param[in]  pSrcB  points to the second input matrix structure
 * @param[out] pDst   points to output matrix structure
 * @return     The function returns either
 * <code>ARM_MATH_SIZE_MISMATCH</code> or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_add_f32(
  const arm_matrix_instance_f32 * pSrcA,
  const arm_matrix_instance_f32 * pSrcB,
        arm_matrix_instance_f32 * pDst);

/**
 * @brief Floating-point matrix subtraction
 * @param[in]  pSrcA  points to the first input matrix structure
 * @param[in]  pSrcB  points to the second input matrix structure
 * @param[out] pDst   points to output matrix structure
 * @return     The function returns either
 * <code>ARM_MATH_SIZE_MISMATCH</code> or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_sub_f32(
  const arm_matrix_instance_f32 * pSrcA,
  const arm_matrix_instance_f32 * pSrcB,
        arm_matrix_instance_f32 * pDst);

/**
 * @brief  Floating-point matrix initialization.
 * @param[in,out] S         points to an instance of the floating-point matrix structure.
 * @param[in]     nRows     number of rows in the matrix.
 * @param[in]     nColumns  number of columns in the matrix.
 * @param[in]     pData     points to the matrix data array.
 */
void arm_mat_init_f32(
        arm_matrix_instance_f32 * S,
        uint16_t nRows,
        uint16_t nColumns,
        float32_t * pData);

/**
 * @brief Floating-point LDL decomposition of Symmetric Positive Semi-Definite Matrix.
 * @param[in]  src   points to the instance of the input floating-point matrix structure.
 * @param[out] l   points to the instance of the output floating-point triangular matrix structure.
 * @param[out] d   points to the instance of the output floating-point diagonal matrix structure.
 * @param[out] p   points to the instance of the output floating-point permutation vector.
 * @return The function returns ARM_MATH_SIZE_MISMATCH, if the dimensions do not match.
 * If the input matrix does not have a decomposition, then the algorithm terminates and returns error status ARM_MATH_DECOMPOSITION_FAILURE.
 * The decomposition is returning a lower triangular matrix.
 */
arm_status arm_mat_ldlt_f32(
  const arm_matrix_instance_f32 * src,
  arm_matrix_instance_f32 * l,
  arm_matrix_instance_f32 * d,
  uint16_t * pp);

/**
 * @brief Floating-point matrix multiplication
 * @param[in]  pSrcA  points to the first input matrix structure
 * @param[in]  pSrcB  points to the second input matrix structure
 * @param[out] pDst   points to output matrix structure
 * @return     The function returns either
 * <code>ARM_MATH_SIZE_MISMATCH</code> or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_mult_f32(
  const arm_matrix_instance_f32 * pSrcA,
  const arm_matrix_instance_f32 * pSrcB,
        arm_matrix_instance_f32 * pDst);

/**
 * @brief Floating-point matrix scaling.
 * @param[in]  pSrc   points to the input matrix
 * @param[in]  scale  scale factor
 * @param[out] pDst   points to the output matrix
 * @return     The function returns either
 * <code>ARM_MATH_SIZE_MISMATCH</code> or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_scale_f32(
  const arm_matrix_instance_f32 * pSrc,
        float32_t scale,
        arm_matrix_instance_f32 * pDst);

/**
 * @brief Floating-point matrix transpose.
 * @param[in]  pSrc  points to the input matrix
 * @param[out] pDst  points to the output matrix
 * @return    The function returns either  <code>ARM_MATH_SIZE_MISMATCH</code>
 * or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_trans_f32(
  const arm_matrix_instance_f32 * pSrc,
        arm_matrix_instance_f32 * pDst);

/**
 * @brief Floating-point matrix and vector multiplication
 * @param[in]  pSrcMat  points to the input matrix structure
 * @param[in]  pVec     points to vector
 * @param[out] pDst     points to output vector
 */
void arm_mat_vec_mult_f32(
  const arm_matrix_instance_f32 *pSrcMat, 
  const float32_t *pVec, 
  float32_t *pDst);

/**
 * @brief         Householder transform of a floating point vector.
 * @param[in]     pSrc        points to the input vector.
 * @param[in]     threshold   norm2 threshold.
 * @param[in]     blockSize   dimension of the vector space.
 * @param[outQ]   pOut        points to the output vector.
 * @return        beta        return the scaling factor beta
 */

float32_t arm_householder_f32(
  const float32_t * pSrc,
  const float32_t threshold,
  uint32_t    blockSize,
  float32_t * pOut);

/**
 * @brief Floating-point matrix inverse.
 * @param[in]  src   points to the instance of the input floating-point matrix structure.
 * @param[out] dst   points to the instance of the output floating-point matrix structure.
 * @return The function returns ARM_MATH_SIZE_MISMATCH, if the dimensions do not match.
 * If the input matrix is singular (does not have an inverse), then the algorithm terminates and returns error status ARM_MATH_SINGULAR.
 */
arm_status arm_mat_inverse_f32(
  const arm_matrix_instance_f32 * src,
  arm_matrix_instance_f32 * dst);

/**
 * @brief Floating-point Cholesky decomposition of Symmetric Positive Definite Matrix.
 * @param[in]  src   points to the instance of the input floating-point matrix structure.
 * @param[out] dst   points to the instance of the output floating-point matrix structure.
 * @return The function returns ARM_MATH_SIZE_MISMATCH, if the dimensions do not match.
 * If the input matrix does not have a decomposition, then the algorithm terminates and returns error status ARM_MATH_DECOMPOSITION_FAILURE.
 * If the matrix is ill conditioned or only semi-definite, then it is better using the LDL^t decomposition.
 * The decomposition is returning a lower triangular matrix.
 */
arm_status arm_mat_cholesky_f32(
  const arm_matrix_instance_f32 * src,
  arm_matrix_instance_f32 * dst);

/**
 * @brief Floating-point, complex, matrix multiplication.
 * @param[in]  pSrcA  points to the first input matrix structure
 * @param[in]  pSrcB  points to the second input matrix structure
 * @param[out] pDst   points to output matrix structure
 * @return     The function returns either
 * <code>ARM_MATH_SIZE_MISMATCH</code> or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_cmplx_mult_f32(
  const arm_matrix_instance_f32 * pSrcA,
  const arm_matrix_instance_f32 * pSrcB,
        arm_matrix_instance_f32 * pDst);

/**
 * @brief Floating-point complex matrix transpose.
 * @param[in]  pSrc  points to the input matrix
 * @param[out] pDst  points to the output matrix
 * @return    The function returns either  <code>ARM_MATH_SIZE_MISMATCH</code>
 * or <code>ARM_MATH_SUCCESS</code> based on the outcome of size checking.
 */
arm_status arm_mat_cmplx_trans_f32(
  const arm_matrix_instance_f32 * pSrc,
  arm_matrix_instance_f32 * pDst);

/**
 * @brief Solve LT . X = A where LT is a lower triangular matrix
 * @param[in]  lt  The lower triangular matrix
 * @param[in]  a  The matrix a
 * @param[out] dst The solution X of LT . X = A
 * @return The function returns ARM_MATH_SINGULAR, if the system can't be solved.
 */
arm_status arm_mat_solve_lower_triangular_f32(
  const arm_matrix_instance_f32 * lt,
  const arm_matrix_instance_f32 * a,
        arm_matrix_instance_f32 * dst);

/**
 * @brief Solve UT . X = A where UT is an upper triangular matrix
 * @param[in]  ut  The upper triangular matrix
 * @param[in]  a  The matrix a
 * @param[out] dst The solution X of UT . X = A
 * @return The function returns ARM_MATH_SINGULAR, if the system can't be solved.
 */
arm_status arm_mat_solve_upper_triangular_f32(
  const arm_matrix_instance_f32 * ut,
  const arm_matrix_instance_f32 * a,
        arm_matrix_instance_f32 * dst);

#ifdef   __cplusplus
}
#endif

#endif