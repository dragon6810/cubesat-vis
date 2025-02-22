#include <Geomag/Geomag.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <miscdefs.h>
#include <VecUtils.h>

#include <Geomag/GeomagLib.h>
#include <CoordinateConversions/CoordinateConversions.h>
#include <Testing/Testing.h>

#include <Geomag/GeomagData.h>

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* INTERNAL DEFINES
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#define WMM_FILENAME ("WMM.COF")

#define ATanH(x)	    (0.5 * log((1 + x) / (1 - x)))

/* Radius of the spherical Earth in meters */
#define EARTH_R ((float)6378100)

/* Angular speed of Earth's rotation about its axis, in rad/s */

#define errorcheck_ret(fresult) if (fresult != FR_OK) return fresult
#define errorcheck_goto(fresult, lbl) if (fresult != FR_OK) goto lbl

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* INTERNAL ROUTINES DEFINITIONS
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

static void Geomag_InitializeModel(Geomag_MagneticModel_t * MagneticModel) {
    MagneticModel->CoefficientFileEndDate = 0;
    MagneticModel->EditionDate = 0;
    strcpy(MagneticModel->ModelName, "");
    MagneticModel->SecularVariationUsed = 0;
    MagneticModel->epoch = 0;
    MagneticModel->nMax = 0;
    MagneticModel->nMaxSecVar = 0;

    for (int i = 0; i < CALCULATE_NUMTERMS(MAX_N_MODE); i++)
    {
        MagneticModel->Main_Field_Coeff_G[i] = 0;
        MagneticModel->Main_Field_Coeff_H[i] = 0;
        MagneticModel->Secular_Var_Coeff_G[i] = 0;
        MagneticModel->Secular_Var_Coeff_H[i] = 0;
    }
}

static FRESULT Geomag_ReadMagModel(char* filename, Geomag_MagneticModel_t *MagneticModel)
{
    int i;

    
    Geomag_InitializeModel(MagneticModel);
    MagneticModel->nMax = MAX_N_MODE;
    MagneticModel->nMaxSecVar = MAX_N_MODE;
    for(i=0; i<CALCULATE_NUMTERMS(MAX_N_MODE)+1; i++)
    {
        MagneticModel->Main_Field_Coeff_G[i] = geomag_wmm_elements[i].coeffg;
        MagneticModel->Main_Field_Coeff_H[i] = geomag_wmm_elements[i].coeffh;
        MagneticModel->Secular_Var_Coeff_G[i] = geomag_wmm_elements[i].slopeg;
        MagneticModel->Secular_Var_Coeff_H[i] = geomag_wmm_elements[i].slopeh;
    }

    return FR_OK;
}

// sets up ellipsoid and geoid to earth's values
static void Geomag_SetDefaults(Geomag_Ellipsoid_t *Ellip, Geomag_Geoid_t *Geoid)
{
    Ellip->a = 6378.137; /*semi-major axis of the ellipsoid in */
    Ellip->b = 6356.7523142; /*semi-minor axis of the ellipsoid in */
    
    // flattening is (a - b) / a
    Ellip->fla = 1 / 298.257223563; /* flattening */
    
    arm_sqrt_f32(1 - (Ellip->b * Ellip->b) / (Ellip->a * Ellip->a), &(Ellip->eps)); /*first eccentricity */
    Ellip->epssq = (Ellip->eps * Ellip->eps); /*first eccentricity squared */
   
    // mean radius
    Ellip->re = 6371.2; /* Earth's radius */

    /* Sets EGM-96 model file parameters */
    Geoid->NumbGeoidCols = 1441; /* 360 degrees of longitude at 15 minute spacing */
    Geoid->NumbGeoidRows = 721; /* 180 degrees of latitude  at 15 minute spacing */
    Geoid->NumbHeaderItems = 6; /* min, max lat, min, max long, lat, long spacing*/
    Geoid->ScaleFactor = 4; /* 4 grid cells per degree at 15 minute spacing  */
    Geoid->NumbGeoidElevs = Geoid->NumbGeoidCols * Geoid->NumbGeoidRows;
    Geoid->Geoid_Initialized = 0; /*  Geoid will be initialized only if this is set to zero */
    Geoid->UseGeoid = MAG_USE_GEOID;
}

static void Geomag_ECEFToGeodetic(Geomag_Ellipsoid_t* Ellip, float32_t x, float32_t y, float32_t z, MAGtype_CoordGeodetic* CoordGeodetic)
{
    const int maxiters = 3;

    int i;

    float d, sx, sy, sz;
    float k, biga, bigb;
    float phi, horlen, h, sinphi, bign;

    sx = x / Ellip->a;
    sy = y / Ellip->a;
    sz = z / Ellip->b;

    // < 1: inside ellipsoid, = 1: on surface of ellipsoid, >1: outside ellipsoid
    d = sx * sx + sy * sy + sz * sz;

    CoordGeodetic->lambda =  atan2f(y, x);

    horlen = sqrtf(x * x + y * y);
    phi = atan2f(z, horlen);
    for(i=0; i<maxiters; i++)
    {
        sinphi = sinf(phi);
        bign = Ellip->a / sqrtf(1.0 - Ellip->epssq * sinphi * sinphi);
        phi = atan2f(z + Ellip->epssq * bign * sinphi, horlen);
    }

    h = horlen / cosf(phi) - bign;

    CoordGeodetic->lambda = RAD2DEG(CoordGeodetic->lambda);
    CoordGeodetic->phi = RAD2DEG(phi);
    CoordGeodetic->HeightAboveEllipsoid = CoordGeodetic->HeightAboveGeoid = h;
}

static float32_t Geomag_TimeToYear(time_t* t)
{
    struct tm date;

    // Shared memory guard because gmtime is not thread-safe
    // SharedMem_BeginAccess();
    struct tm *temp = gmtime(t);
    memcpy(&date, temp, sizeof(struct tm));
    // SharedMem_EndAccess();

    int year = date.tm_year + 1900;

    // TODO: Check approximation
    float32_t days = (year % 4 == 0) ? 366.0 : 365.0;

    return year + (date.tm_yday)/(days);
}

void Geomag_RunTests(const char *filename)
{
    int i;

    float year, altitude, latitude, longitude, testp;
    time_t time;

    Geomag_MagneticModel_t MagneticModel, TimedMagneticModel;
    Geomag_Ellipsoid_t Ellip;
    Geomag_CoordSpherical_t CoordSpherical;
    Geomag_Geoid_t Geoid;
    Vec3D_t cart, vec, NED, ECEF, ECI;
    char testname[12];

    Geomag_InitializeModel(&MagneticModel);
    Geomag_ReadMagModel(WMM_FILENAME, &MagneticModel);
    Geomag_InitializeModel(&TimedMagneticModel);
    Geomag_SetDefaults(&Ellip, &Geoid);

    strcpy(testname, "WMM Test 00");

    i = 0;
    for(i=0; i<sizeof(geomag_wmm_test_elements) / sizeof(wmm_test_element_t); i++)
    {
        year = geomag_wmm_test_elements[i].time;
        altitude = geomag_wmm_test_elements[i].alt;
        latitude = geomag_wmm_test_elements[i].lat;
        longitude = geomag_wmm_test_elements[i].lon;
        NED.Vec[0] = geomag_wmm_test_elements[i].x;
        NED.Vec[1] = geomag_wmm_test_elements[i].y;
        NED.Vec[2] = geomag_wmm_test_elements[i].z;

        time = (time_t)((year - 1970) * 365.25 * 24 * 3600);

        CoordSpherical.r = Ellip.re + altitude;
        CoordSpherical.phig = DEG2RAD(latitude);
        CoordSpherical.lambda = DEG2RAD(longitude);

        Coord_NEDToECEF(&NED, &ECEF, latitude, longitude);
        Coord_ECEFToECI(time, &ECEF, &ECI);

        CoordSpherical.lambda;
        cart.X = CoordSpherical.r * cosf(CoordSpherical.lambda) * cosf(CoordSpherical.phig);
        cart.Y = CoordSpherical.r * sinf(CoordSpherical.lambda) * cosf(CoordSpherical.phig);
        cart.Z = CoordSpherical.r * sinf(CoordSpherical.phig);

        Coord_ECEFToECI(time, &cart, &cart);
        Geomag_GetMagEquatorial(&time, &cart, &vec);

        assert(i < 100);
        testname[9] = (i / 10) + '0';
        testname[10] = (i % 10) + '0';

        printf("%s:\n   { %f, %f, %f }\n    { %f, %f, %f }\n", testname, ECI.X, ECI.Y, ECI.Z, vec.X, vec.Y, vec.Z);

        Testing_TestVector(ECI, cart, vec, 512.0, testname);
    }
}

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* EXTERNAL ROUTINES DEFINITIONS
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
    Take in ECI coordinates and spits out ECI coordinates
*/
FRESULT Geomag_GetMagEquatorial(time_t* t, const Vec3D_t* SatECI, Vec3D_t* MagEquatorial)
{
    Vec3D_t SatECEF;
    Vec3D_t MagECEF;
    float32_t theta, phi;
    Geomag_MagneticModel_t MagneticModel;
    Geomag_MagneticModel_t TimedMagneticModel;
    MAGtype_Ellipsoid Ellip;
    Geomag_Geoid_t Geoid;
    MAGtype_CoordGeodetic CoordGeodetic;
    MAGtype_CoordSpherical CoordSpherical;
    MAGtype_GeoMagneticElements GeomagElements;
    Vec3D_t MagNED;
    MAGtype_GeoMagneticElements GeomagErrors;
    MAGtype_Date date;
    
    Coord_ECIToECEF(*t, SatECI, &SatECEF);

    // theta is latitude, phi is longitude.
    arm_atan2_f32(SatECEF.Z, sqrtf(SatECEF.X*SatECEF.X + SatECEF.Y*SatECEF.Y), &theta);
    arm_atan2_f32(SatECEF.Y, SatECEF.X, &phi);
    theta = RAD2DEG(theta);
    phi = RAD2DEG(phi);
    
    Geomag_SetDefaults(&Ellip, &Geoid);
    Geoid.Geoid_Initialized = pdTRUE;

    Geomag_ECEFToGeodetic(&Ellip, SatECEF.X, SatECEF.Y, SatECEF.Z, &CoordGeodetic);

    date.DecimalYear = Geomag_TimeToYear(t);
    calculateMagneticField(&CoordGeodetic, &date, &GeomagElements, &GeomagErrors);

    MagNED.X = GeomagElements.X;
    MagNED.Y = GeomagElements.Y;
    MagNED.Z = GeomagElements.Z;
    Coord_NEDToECEF(&MagNED, &MagECEF, theta, phi);
    Coord_ECEFToECI(*t, &MagECEF, MagEquatorial);

    return FR_OK;
}
