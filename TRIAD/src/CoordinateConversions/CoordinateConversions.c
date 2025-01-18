#include <CoordinateConversions/CoordinateConversions.h>

#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>

#include <VecUtils.h>

#include <Testing/Testing.h>

void Coord_NEDToECEF(Vec3D_t* NED, Vec3D_t *ECEF, float32_t theta, float32_t phi_p)
{
    float lambda, phi;
    float sinlambda, coslambda;
    float sinphi, cosphi;
    Vec3D_t spherical;

    lambda = DEG2RAD(phi_p);
    phi = DEG2RAD(theta);

    sinlambda = sinf(lambda);
    coslambda = cosf(lambda);
    sinphi = sinf(phi);
    cosphi = cosf(phi);

    ECEF->X =
        NED->X * -sinphi *  coslambda + 
        NED->Y *  1.0    * -sinlambda + 
        NED->Z * -cosphi *  coslambda ;

    ECEF->Y =
        NED->X * -sinphi *  sinlambda + 
        NED->Y *  1.0    *  coslambda + 
        NED->Z * -cosphi *  sinlambda ;

    ECEF->Z =
        NED->X *  cosphi *  1.0       + 
        NED->Y *  0.0    *  0.0       + 
        NED->Z * -sinphi *  1.0       ;
}

void Coord_ECEFToECI(time_t t, const Vec3D_t* ECEF, Vec3D_t* ECI)
{
    assert(ECEF);
    assert(ECI);

    // The exact negative of Geomac_ECIToECEF
    float32_t rotAngle = EARTH_ROT_SPEED * (t - EQUINOX_TIME);
    
    Vec_RotateSpher(ECEF, 0, rotAngle, ECI);
}

void Coord_ECIToECEF(time_t t, const Vec3D_t* ECI, Vec3D_t* ECEF)
{
    assert(ECI);
    assert(ECEF);

    // Angle between prime meridian and vernal equinox
    // Distance we'd need to cover to align our prime meridian with vernal equinox.
    // rotation = speed * time so
    // rotAngle = EARTH_ROT_SPEED * EQUINOX_TIME - EARTH_ROT_SPEED * t so
    // rotAngle = EARTH_ROT_SPEED * (EQUINOX_TIME - t);
    float32_t rotAngle = EARTH_ROT_SPEED * (EQUINOX_TIME - t);
    
    Vec_RotateSpher(ECI, 0, rotAngle, ECEF);
}

static void Coord_TestNEDToECEF(void)
{
    Vec3D_t NED, ECEF, expect;
    float theta, phi;

    // Nominal cases

    NED.X = 1000.0; NED.Y = 2000.0; NED.Z = -500.0;
    theta = 0;
    phi = 0;
    expect.X = 500.0; expect.Y = 2000.0; expect.Z = 1000.0;
    Coord_NEDToECEF(&NED, &ECEF, theta, phi);
    Testing_TestVector(expect, NED, ECEF, 0.01, "Coord_NEDToECEF Nominal Test 1");

    NED.X = 1000.0; NED.Y = 2000.0; NED.Z = -500;
    theta = 0;
    phi = 180;
    expect.X = -500.0; expect.Y = -2000.0; expect.Z = 1000.0;
    Coord_NEDToECEF(&NED, &ECEF, theta, phi);
    Testing_TestVector(expect, NED, ECEF, 0.01, "Coord_NEDToECEF Nominal Test 2");
}

void Coord_TestConversions(void)
{
    Coord_TestNEDToECEF();
}