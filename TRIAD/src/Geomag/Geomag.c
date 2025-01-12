#include <Geomag/Geomag.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <miscdefs.h>
#include <VecUtils.h>

#include <Geomag/GeomagLib.h>

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* INTERNAL DEFINES
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#define WMM_FILENAME ("WMM.COF")

#define ATanH(x)	    (0.5 * log((1 + x) / (1 - x)))

/* UTC time at which 0deg N and 0deg E was along the X axis in equatorial coordinates.
   This happens at the time of a vernal (spring) equinox.
   (Actually this also happens every day whenever RA 0 crosses the meridian at 0deg N 0deg E,
   but an equinox is a more convenient reference point.) */
#define EQUINOX_TIME ((time_t) 1695451800)
/* Radius of the spherical Earth in meters */
#define EARTH_R ((float)6378100)

/* Angular speed of Earth's rotation about its axis, in rad/s */
#define EARTH_ROT_SPEED ((float)7.2921159E-5)

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

    for (int i = 0; i < CALCULATE_NUMTERMS(MAX_N_MODE); i++) {
        MagneticModel->Main_Field_Coeff_G[i] = 0;
        MagneticModel->Main_Field_Coeff_H[i] = 0;
        MagneticModel->Secular_Var_Coeff_G[i] = 0;
        MagneticModel->Secular_Var_Coeff_H[i] = 0;
    }
}

static FRESULT Geomag_ReadMagModel(char* filename, Geomag_MagneticModel_t *MagneticModel)
{
    // Type to read each line into
    typedef struct {
        int n, m;
        float32_t gnm, hnm, dgnm, dhnm;
    } WMM_Line;

    FILE *ptr;

    // Setup variable
    Geomag_InitializeModel(MagneticModel);
    WMM_Line line;
    int index;
    int count;
    int bogus;

    ptr = fopen(WMM_FILENAME, "r");
    if(!ptr)
    {
        printf("can't open file \"%s\".\n", WMM_FILENAME);
        return FR_NO_FILE;
    }

    fscanf(ptr, " %f WMM-%d %d/%d/%d\n", &MagneticModel->epoch, &bogus, &bogus, &bogus, &bogus);

    MagneticModel->Main_Field_Coeff_H[0] = 0.0;
    MagneticModel->Main_Field_Coeff_G[0] = 0.0;
    MagneticModel->Secular_Var_Coeff_H[0] = 0.0;
    MagneticModel->Secular_Var_Coeff_G[0] = 0.0;

    count = 0;
    while(fscanf(ptr, " %d %d %f %f %f %f", &line.n, &line.m, &line.gnm, &line.hnm, &line.dgnm, &line.dhnm) == 6)
    {
        count++;
        assert(count < CALCULATE_NUMTERMS(MAX_N_MODE)+1);

        MagneticModel->Main_Field_Coeff_G[count] = line.gnm;
        MagneticModel->Main_Field_Coeff_H[count] = line.hnm;
        MagneticModel->Secular_Var_Coeff_G[count] = line.dgnm;
        MagneticModel->Secular_Var_Coeff_H[count] = line.dhnm;
    }

    fclose(ptr);

    MagneticModel->nMax = MAX_N_MODE;
    MagneticModel->nMaxSecVar = MAX_N_MODE;

    return FR_OK;
}

// sets up ellipsoid and geoid to earth's values
static void Geomag_SetDefaults(Geomag_Ellipsoid_t *Ellip, Geomag_Geoid_t *Geoid) {
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

static void Geomag_SphericalToCartesian(Geomag_CoordSpherical_t* Sphe, Vec3D_t* Cart)
{
    Cart->X = Sphe->r * cosf(Sphe->lambda) * cosf(Sphe->phig);
    Cart->Y = Sphe->r * sinf(Sphe->lambda) * cosf(Sphe->phig);
    Cart->Z = Sphe->r * sinf(Sphe->phig  );
}

static void Geomag_CartesianToGeodetic(Geomag_Ellipsoid_t* Ellip, float32_t x, float32_t y, float32_t z, MAGtype_CoordGeodetic* CoordGeodetic)
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

static float32_t Geomag_TimeToYear(time_t* t) {
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

static void Geomag_NEDToECEF(MAGtype_GeoMagneticElements* MagHorizontal, Vec3D_t *MagEquatorial, float32_t theta, float32_t phi_p)
{
    float lambda, phi;
    Vec3D_t spherical;

    lambda = DEG2RAD(phi_p);
    phi = DEG2RAD(theta);

    MagEquatorial->X =
        MagHorizontal->X * -sinf(phi) *  cosf(lambda) + 
        MagHorizontal->Y *  1.0       * -sinf(lambda) + 
        MagHorizontal->Z * -cosf(phi) *  cosf(lambda);

    MagEquatorial->Y =
        MagHorizontal->X * -sinf(phi) *  sinf(lambda) + 
        MagHorizontal->Y *  1.0       *  cosf(lambda) + 
        MagHorizontal->Z * -cosf(phi) *  sinf(lambda);

    MagEquatorial->Z =
        MagHorizontal->X *  cosf(phi) *  1.0          + 
        MagHorizontal->Y *  0.0       *  0.0          + 
        MagHorizontal->Z * -sinf(phi) *  1.0;
}

static void Geomag_TestHorizontalToEquatorial(void)
{
    Vec3D_t MagEquatorial;
    MAGtype_GeoMagneticElements MagHorizontal;

    MagHorizontal.X = MagHorizontal.Y = MagHorizontal.Z = 0;
    MagHorizontal.X = 20;
    MagHorizontal.Y = 0;
    MagHorizontal.Z = 1;

    Geomag_NEDToECEF(&MagHorizontal, &MagEquatorial, 0, 90);
    printf("Mag Equatorial: %f, %f, %f.\n", MagEquatorial.X, MagEquatorial.Y, MagEquatorial.Z);
}

// math mostly from https://journals.pan.pl/Content/98324/PDF/art05.pdf
static void Geomag_TestCartesianToGeodetic(void)
{
    MAGtype_Ellipsoid Ellip;
    Geomag_Geoid_t Geoid;
    Vec3D_t cartesian;
    MAGtype_CoordGeodetic geodetic;

    Geomag_SetDefaults(&Ellip, &Geoid);
    Geoid.Geoid_Initialized = pdTRUE;

    cartesian.X = 7000;
    cartesian.Y = -7000;
    cartesian.Z = 7000;

    Geomag_CartesianToGeodetic(&Ellip, cartesian.X, cartesian.Y, cartesian.Z, &geodetic);
    printf("Coord Cartesian: %f, %f, %f.\n", cartesian.X, cartesian.Y, cartesian.Z);
    printf("Coord Geodetic: %f, %f, %f.\n", geodetic.lambda, geodetic.phi, geodetic.HeightAboveEllipsoid);
}

void Geomag_RunTests(const char *filename)
{
    int i;

    FILE *ptr;
    float year, altitude, latitude, longitude, testp, x, y, z;
    time_t time;
    float p;
    float delta;

    ptr = fopen(filename, "r");
    if (!ptr)
    {
        printf("couldn't open test file \"%s\".\n", filename);
        abort();
    }

    // skip the spec at the top by skipping 18 lines
    char buffer[256];
    for (i=0; i<18; ++i)
        fgets(buffer, sizeof(buffer), ptr);

    Geomag_MagneticModel_t MagneticModel, TimedMagneticModel;
    Geomag_Ellipsoid_t Ellip;
    Geomag_CoordSpherical_t CoordSpherical;
    Geomag_Geoid_t Geoid;
    Vec3D_t cart, vec;

    Geomag_InitializeModel(&MagneticModel);
    Geomag_ReadMagModel(WMM_FILENAME, &MagneticModel);
    Geomag_InitializeModel(&TimedMagneticModel);
    Geomag_SetDefaults(&Ellip, &Geoid);

    i = 0;
    while (fscanf(ptr, "%f %f %f %f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f\n", 
                  &year, &altitude, &latitude, &longitude, &x, &y, &z, &testp) == 8) 
    {
        time = (time_t)((year - 1970) * 365.25 * 24 * 3600);

        CoordSpherical.r = Ellip.re + altitude;
        CoordSpherical.phig = DEG2RAD(latitude);
        CoordSpherical.lambda = DEG2RAD(longitude) + EARTH_ROT_SPEED * (time - EQUINOX_TIME); // relative to vernal equinox

        Geomag_SphericalToCartesian(&CoordSpherical, &cart);

        cart.X = latitude;
        cart.Y = longitude;
        cart.Z = altitude;
        Geomag_GetMagEquatorial(&time, &cart, &vec);
        p = sqrtf(vec.X * vec.X + vec.Y * vec.Y + vec.Z * vec.Z);
        testp = sqrtf(x * x + y * y + z * z);

        printf("test %d:\n", i++);
        printf("  Input: Year=%.2f, Altitude=%.2f km, Latitude=%.2f deg, Longitude=%.2f deg\n",
               year, altitude, latitude, longitude);
        printf("  Expected Potential: %.6f\n", testp);
        printf("  Computed Potential: %.6f\n", p);

        delta = fabs(p - testp);
        printf("  Difference: %.6f\n\n", delta);
    }

    Geomag_TestHorizontalToEquatorial();
    Geomag_TestCartesianToGeodetic();

    fclose(ptr);
}

void Geomag_ECEFToECI(time_t t, const Vec3D_t* ECEF, Vec3D_t* ECI)
{
    assert(ECEF);
    assert(ECI);

    // The exact negative of Geomac_ECIToECEF
    float32_t rotAngle = EARTH_ROT_SPEED * (t - EQUINOX_TIME);
    
    Vec_RotateSpher(ECEF, 0, rotAngle, ECI);
}

void Geomag_ECIToECEF(time_t t, const Vec3D_t* ECI, Vec3D_t* ECEF)
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
    MAGtype_GeoMagneticElements GeomagErrors;
    MAGtype_Date date;
    
    Geomag_ECIToECEF(*t, SatECI, &SatECEF);
    
    // theta is latitude, phi is longitude.
    arm_atan2_f32(sqrtf(SatECEF.X*SatECEF.X + SatECEF.Y*SatECEF.Y), SatECEF.Z, &theta);
    arm_atan2_f32(SatECEF.Y, SatECEF.X, &phi);
    theta = RAD2DEG(theta);
    phi = RAD2DEG(phi);
    printf("lat, long: (%f, %f).\n", theta, phi);
    
    Geomag_SetDefaults(&Ellip, &Geoid);
    Geoid.Geoid_Initialized = pdTRUE;

    Geomag_CartesianToGeodetic(&Ellip, SatECEF.X, SatECEF.Y, SatECEF.Z, &CoordGeodetic);

    date.DecimalYear = Geomag_TimeToYear(t);

    calculateMagneticField(&CoordGeodetic, &date, &GeomagElements, &GeomagErrors);

    Geomag_NEDToECEF(&GeomagElements, &MagECEF, theta, phi);
    Geomag_ECEFToECI(*t, &MagECEF, MagEquatorial);

    return FR_OK;
}
