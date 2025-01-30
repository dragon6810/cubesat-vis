#include <arm_math.h>

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

arm_status arm_mat_sub_f32(
  const arm_matrix_instance_f32 * pSrcA,
  const arm_matrix_instance_f32 * pSrcB,
        arm_matrix_instance_f32 * pDst)
{
    int i;

    assert(pSrcA);
    assert(pSrcB);
    assert(pDst);
    assert(pSrcA->pData);
    assert(pSrcB->pData);
    assert(pDst->pData);

    if(pSrcA->numCols != pSrcB->numCols)
        return ARM_MATH_SIZE_MISMATCH;

    if(pSrcA->numRows != pSrcB->numRows)
        return ARM_MATH_SIZE_MISMATCH;

    if(pSrcA->numCols != pDst->numCols)
        return ARM_MATH_SIZE_MISMATCH;

    if(pSrcA->numRows != pDst->numRows)
        return ARM_MATH_SIZE_MISMATCH;

    for(i=0; i<pSrcA->numCols*pSrcA->numRows; i++)
        pDst->pData[i] = pSrcA->pData[i] - pSrcB->pData[i];

    return ARM_MATH_SUCCESS;
}

void arm_mat_init_f32(
        arm_matrix_instance_f32 * S,
        uint16_t nRows,
        uint16_t nColumns,
        float32_t * pData)
{
    assert(S);

    S->numRows = nRows;
    S->numCols = nColumns;
    S->pData = pData;
}

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
        arm_matrix_instance_f32 * pDst)
{
    int i, j, k;

    int cura, curb, curc;

    assert(pSrcA);
    assert(pSrcB);
    assert(pDst);

    assert(pSrcA->numCols == pSrcB->numRows);
    assert(pSrcA->numRows == pDst->numRows);
    assert(pSrcB->numCols == pDst->numCols);

    for(i=0; i<pSrcA->numRows; i++)
    {
        for(j=0; j<pSrcB->numCols; j++)
        {
            curc = i*pDst->numCols+j;
            pDst->pData[curc] = 0;
            for(k=0; k<pSrcB->numRows; k++)
            {
                cura = i*pSrcA->numCols+k;
                curb = k*pSrcB->numCols+j;

                pDst->pData[curc] += pSrcA->pData[cura] * pSrcB->pData[curb];
            }
        }
    }

#if 0
    printf("\n");
    printf("[ %f %f %f ]   [ %f %f %f ]   [ %f %f %f ]\n", 
        pSrcA->pData[0], pSrcA->pData[1], pSrcA->pData[2], 
        pSrcB->pData[0], pSrcB->pData[1], pSrcB->pData[2],
         pDst->pData[0],  pDst->pData[1],  pDst->pData[2]);
    printf("[ %f %f %f ] X [ %f %f %f ] = [ %f %f %f ]\n", 
        pSrcA->pData[3], pSrcA->pData[4], pSrcA->pData[5], 
        pSrcB->pData[3], pSrcB->pData[4], pSrcB->pData[5],
         pDst->pData[3],  pDst->pData[4],  pDst->pData[5]);
    printf("[ %f %f %f ]   [ %f %f %f ]   [ %f %f %f ]\n", 
        pSrcA->pData[6], pSrcA->pData[7], pSrcA->pData[8], 
        pSrcB->pData[6], pSrcB->pData[7], pSrcB->pData[8],
         pDst->pData[6],  pDst->pData[7],  pDst->pData[8]);
#endif

    return ARM_MATH_SUCCESS;
}

arm_status arm_mat_scale_f32(
  const arm_matrix_instance_f32 * pSrc,
        float32_t scale,
        arm_matrix_instance_f32 * pDst)
{
    int i;

    assert(pSrc);
    assert(pDst);
    assert(pSrc->pData);
    assert(pDst->pData);

    if(pSrc->numCols != pDst->numCols)
        return ARM_MATH_SIZE_MISMATCH;
    if(pSrc->numRows != pDst->numRows)
        return ARM_MATH_SIZE_MISMATCH;

    for(i=0; i<pSrc->numCols*pSrc->numRows; i++)
        pDst->pData[i] = pSrc->pData[i] * scale;

    return ARM_MATH_SUCCESS;
}

void arm_mat_vec_mult_f32(
  const arm_matrix_instance_f32 *pSrcMat, 
  const float32_t *pVec, 
  float32_t *pDst)
{
    int i, j;

    float32_t *mat;

    assert(pSrcMat);
    assert(pSrcMat->pData);
    assert(pVec);

    float32_t v[pSrcMat->numRows];

    for(i=0, mat=pSrcMat->pData; i<pSrcMat->numRows; i++, mat += pSrcMat->numCols)
    {
        v[i] = 0;
        for(j=0; j<pSrcMat->numCols; j++)
            v[i] += pVec[j] * mat[j];
    }

    if(pDst)
        memcpy(pDst, v, pSrcMat->numRows * sizeof(float32_t));
    else
        memcpy(pVec, v, pSrcMat->numRows * sizeof(float32_t));
}