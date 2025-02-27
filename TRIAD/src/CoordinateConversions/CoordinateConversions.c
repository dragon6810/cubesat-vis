#include <CoordinateConversions/CoordinateConversions.h>

#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>

#include <VecUtils.h>

#include <Testing/Testing.h>

void Coord_NEDToECEF(Vec3D_t* NED, Vec3D_t *ECEF, float32_t lat, float32_t lon)
{
    float lambda, phi;
    float sinlambda, coslambda;
    float sinphi, cosphi;
    Vec3D_t spherical;

    lambda = DEG2RAD(lon);
    phi = DEG2RAD(lat);

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
    float lat, lon;

    // Nominal cases

    NED.X = 1000.0; NED.Y = 2000.0; NED.Z = -500.0;
    lat = 0;
    lon = 0;
    expect.X = 500.0; expect.Y = 2000.0; expect.Z = 1000.0;
    Coord_NEDToECEF(&NED, &ECEF, lat, lon);
    Testing_TestVector(expect, NED, ECEF, 0.01, "Coord_NEDToECEF Nominal Test 1");

    NED.X = 1000.0; NED.Y = 2000.0; NED.Z = -500;
    lat = 0;
    lon = 180;
    expect.X = -500.0; expect.Y = -2000.0; expect.Z = 1000.0;
    Coord_NEDToECEF(&NED, &ECEF, lat, lon);
    Testing_TestVector(expect, NED, ECEF, 0.01, "Coord_NEDToECEF Nominal Test 2");

    NED.X = 1000.0; NED.Y = 1000.0; NED.Z = 0;
    lat = 45;
    lon = 90;
    expect.X = -1000.0; expect.Y = -707.10678; expect.Z = 707.10678;
    Coord_NEDToECEF(&NED, &ECEF, lat, lon);
    Testing_TestVector(expect, NED, ECEF, 0.01, "Coord_NEDToECEF Nominal Test 3");
}

void Coord_TestConversions(void)
{
    Coord_TestNEDToECEF();
}