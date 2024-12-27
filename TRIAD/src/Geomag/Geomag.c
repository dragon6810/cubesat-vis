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
    float32_t b,r,e,f,p,q,d,v,g,t,zlong,rlat;

    /*
    *   1.0 compute semi-minor axis and set sign to that of z in order
    *       to get sign of Phi correct
    */

    if (z < 0.0)
        b = -Ellip->b;
    else
        b = Ellip->b;

    /*
    *   2.0 compute intermediate values for latitude
    */
    arm_sqrt_f32(x*x + y*y, &r);
    e = ( b*z - (Ellip->a*Ellip->a - b*b) ) / ( Ellip->a*r );
    f = ( b*z + (Ellip->a*Ellip->a - b*b) ) / ( Ellip->a*r );
    
    /*
    *   3.0 find solution to:
    *       t^4 + 2*E*t^3 + 2*F*t - 1 = 0
    */
    p= (4.0 / 3.0) * (e*f + 1.0);
    q= 2.0 * (e*e - f*f);
    d= p*p*p + q*q;

    if( d >= 0.0 )
    {
        arm_sqrt_f32(d, &d);
        v= pow( (d - q), (1.0 / 3.0) ) - pow( (d+ q), (1.0 / 3.0) );
    }
    else
        {
        float32_t sqrt_p;
        arm_sqrt_f32(-p, &sqrt_p);
        v= 2.0 * sqrt_p
            * arm_cos_f32( acos( q/(p * sqrt_p) ) / 3.0 );
        }
    /*
    *   4.0 improve v
    *       NOTE: not really necessary unless point is near pole
    */
    if( v*v < fabs(p) ) {
            v= -(v*v*v + 2.0*q) / (3.0*p);
    }
    arm_sqrt_f32(e*e + v, &g);
    g = (g + e) / 2.0;
    arm_sqrt_f32(g*g  + (f - v*g)/(2.0*g - e), &t);
    t = t - g;

    arm_atan2_f32((Ellip->a * (1.0 - t * t)), (2.0 * b * t), &rlat); // ORIGINAL: rlat =atan( (Ellip->a*(1.0 - t*t)) / (2.0*b*t) );
    CoordGeodetic->phi = RAD2DEG(rlat);

    /*
    *   5.0 compute height above ellipsoid
    */
    CoordGeodetic->HeightAboveEllipsoid = (r - Ellip->a*t) * arm_cos_f32(rlat) + (z - b) * arm_sin_f32(rlat);

    /*
    *   6.0 compute longitude east of Greenwich
    */
    arm_atan2_f32(y, x, &zlong);
    if( zlong < 0.0 )
            zlong= zlong + 2*PI;

    CoordGeodetic->lambda = RAD2DEG(zlong);
    while(CoordGeodetic->lambda > 180)
    {
        CoordGeodetic->lambda-=360;
    }
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

static void Geomag_HorizontalToEquatorial(MAGtype_GeoMagneticElements* MagHorizontal, Vec3D_t *MagEquatorial, float32_t theta, float32_t phi_p)
{
    float lambda, phi;
    Vec3D_t spherical;

    lambda = DEG2RAD(phi_p);
    phi = DEG2RAD(theta);

    MagEquatorial->X =
        MagHorizontal->X * -sinf(phi) * cosf(lambda) + 
        MagHorizontal->Y *  1.0       * sinf(lambda) + 
        MagHorizontal->Z * -cosf(phi) * cosf(lambda);

    MagEquatorial->Y =
        MagHorizontal->X * -sinf(phi) * sinf(lambda) + 
        MagHorizontal->Y *  1.0       * cosf(lambda) + 
        MagHorizontal->Z * -cosf(phi) * sinf(lambda);

    MagEquatorial->Z =
        MagHorizontal->X *  sinf(phi) * 1.0          + 
        MagHorizontal->Y *  0.0       * 0.0          + 
        MagHorizontal->Z * -sinf(phi) * 1.0;
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

    fclose(ptr);
}

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* EXTERNAL ROUTINES DEFINITIONS
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

FRESULT Geomag_GetMagEquatorial(time_t* t, const Vec3D_t* SatEquatorial, Vec3D_t* MagEquatorial)
{
    Vec3D_t SatLocal;
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

    // Angle between prime meridian and vernal equinox
    // Distance we'd need to cover to align our prime meridian with vernal equinox.
    // rotation = speed * time so
    // rotAngle = EARTH_ROT_SPEED * EQUINOX_TIME - EARTH_ROT_SPEED * t so
    // rotAngle = EARTH_ROT_SPEED * (EQUINOX_TIME - t);
    float32_t rotAngle = EARTH_ROT_SPEED * (EQUINOX_TIME - (*t));
    // SatLocal is relative to prime meridian
    Vec_RotateSpher(SatEquatorial, 0, rotAngle, &SatLocal);
    
    // Lattitude and Longitude angles for rotation at the end
    // theta is longitude, phi is latitude.
    arm_atan2_f32(sqrtf(SatLocal.X*SatLocal.X + SatLocal.Y*SatLocal.Y), SatLocal.Z, &theta);
    arm_atan2_f32(SatLocal.Y, SatLocal.X, &phi);
    
    Geomag_SetDefaults(&Ellip, &Geoid);
    Geoid.Geoid_Initialized = pdTRUE;

    Geomag_CartesianToGeodetic(&Ellip, SatLocal.X, SatLocal.Y, SatLocal.Z, &CoordGeodetic);
    // Geomag_GeodeticToSpherical(&Ellip, &CoordGeodetic, &CoordSpherical);

    date.DecimalYear = Geomag_TimeToYear(t);

    CoordGeodetic.lambda = RAD2DEG(CoordGeodetic.lambda);
    CoordGeodetic.phi = RAD2DEG(CoordGeodetic.phi);

    CoordGeodetic.phi = -CoordGeodetic.phi + 90;

    CoordGeodetic.phi = SatEquatorial->X;
    CoordGeodetic.lambda = SatEquatorial->Y;
    CoordGeodetic.HeightAboveEllipsoid = SatEquatorial->Z;

    while(CoordGeodetic.lambda < -360)
        CoordGeodetic.lambda += 360;
    while(CoordGeodetic.lambda > 360)
        CoordGeodetic.lambda -= 360;

    while(CoordGeodetic.phi < -360)
        CoordGeodetic.phi += 360;
    while(CoordGeodetic.phi > 360)
        CoordGeodetic.phi -= 360;

    calculateMagneticField(&CoordGeodetic, &date, &GeomagElements, &GeomagErrors);
    MagEquatorial->X = GeomagElements.X;
    MagEquatorial->Y = GeomagElements.Y;
    MagEquatorial->Z = GeomagElements.Z;

    theta = CoordGeodetic.phi;
    phi = CoordGeodetic.lambda;

    Geomag_HorizontalToEquatorial(&GeomagElements, MagEquatorial, theta, phi);

    return FR_OK;
}
