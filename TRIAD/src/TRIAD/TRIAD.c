#include <TRIAD/TRIAD.h>

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

static void TRIAD_MatTrans3x3(arm_matrix_instance_f32* mat)
{
    int i, j;

    float old[9];

    assert(mat->numCols == 3);
    assert(mat->numRows == 3);
    assert(mat->pData);

    memcpy(old, mat->pData, 9 * sizeof(float));

    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
            mat->pData[i * 3 +j] = old[j * 3 + i];
    }
}

static Vec3D_t TRIAD_CrossVec3(Vec3D_t a, Vec3D_t b)
{
    Vec3D_t v;

    v.Vec[0] = a.Vec[1] * b.Vec[2] - a.Vec[2] * b.Vec[1];
    v.Vec[1] = a.Vec[2] * b.Vec[0] - a.Vec[0] * b.Vec[2];
    v.Vec[2] = a.Vec[0] * b.Vec[1] - a.Vec[1] * b.Vec[0];

    return v;
}

static Vec3D_t TRIAD_NormalizeVec3(Vec3D_t v)
{
    int i;

    Vec3D_t newv;
    float len;

    for(i=0, len=0.0; i<3; i++)
        len += v.Vec[i] * v.Vec[i];
    len = sqrtf(len);

    for(i=0; i<3; i++)
        newv.Vec[i] = v.Vec[i] / len;

    return newv;
}

int TRIAD_Compute(Vec3D_t* r1, Vec3D_t* r2, Vec3D_t* R1, Vec3D_t* R2, arm_matrix_instance_f32* mat)
{
    int i;

    Vec3D_t s, S, m, M, sxm, SxM;
    arm_matrix_instance_f32 a, b, att;
    float32_t adata[9], bdata[9], attdata[9];

    assert(mat);
    assert(mat->numRows == 3);
    assert(mat->numCols == 3);
    assert(mat->pData);

    s = TRIAD_NormalizeVec3(*r1);
    S = TRIAD_NormalizeVec3(*R1);

    m = TRIAD_NormalizeVec3(TRIAD_CrossVec3(*r1, *r2));
    M = TRIAD_NormalizeVec3(TRIAD_CrossVec3(*R1, *R2));

    sxm = TRIAD_CrossVec3(s, m);
    SxM = TRIAD_CrossVec3(S, M);

    arm_mat_init_f32(&a, 3, 3, adata);
    arm_mat_init_f32(&b, 3, 3, bdata);
    arm_mat_init_f32(&att, 3, 3, attdata);

    adata[0] = S.X; adata[1] = M.X; adata[2] = SxM.X;
    adata[3] = S.Y; adata[4] = M.Y; adata[5] = SxM.Y;
    adata[6] = S.Z; adata[7] = M.Z; adata[8] = SxM.Z;

    bdata[0] = s.X;   bdata[1] = s.Y;   bdata[2] = s.Z;
    bdata[3] = m.X;   bdata[4] = m.Y;   bdata[5] = m.Z;
    bdata[6] = sxm.X; bdata[7] = sxm.Y; bdata[8] = sxm.Z;

    arm_mat_mult_f32(&a, &b, &att);

    memcpy(mat->pData, attdata, 9 * sizeof(float32_t));

    return 1;
}

int TRIAD_CalculateOmega(arm_matrix_instance_f32* atta, arm_matrix_instance_f32* attb, float interval, Vec3D_t* omega)
{
    arm_matrix_instance_f32 diff;
    arm_matrix_instance_f32 attt;
    arm_matrix_instance_f32 omegamat;
    float32_t diffdata[3 * 3];
    float32_t atttdata[3 * 3];
    float32_t omegamatdata[3 * 3];

    assert(atta);
    assert(attb);
    assert(omega);
    assert(atta->numCols == 3);
    assert(atta->numRows == 3);
    assert(attb->numCols == 3);
    assert(attb->numRows == 3);
    assert(atta->pData);
    assert(attb->pData);

    arm_mat_init_f32(&diff, 3, 3, diffdata);
    arm_mat_init_f32(&attt, 3, 3, atttdata);
    arm_mat_init_f32(&omegamat, 3, 3, omegamatdata);
    
    arm_mat_sub_f32(atta, attb, &diff);
    arm_mat_scale_f32(&diff, 1.0 / interval, &diff);
    
    memcpy(attt.pData, atta->pData, 3 * 3 * sizeof(float));
    TRIAD_MatTrans3x3(&attt);

    arm_mat_mult_f32(&diff, &attt, &omegamat);

    omega->X = -omegamat.pData[7];
    omega->Y = -omegamat.pData[2];
    omega->Z = -omegamat.pData[3];

    return 1;
}