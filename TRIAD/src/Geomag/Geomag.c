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

static void Geomag_GeodeticToSpherical(Geomag_Ellipsoid_t* Ellip, MAGtype_CoordGeodetic* CoordGeodetic, MAGtype_CoordSpherical *CoordSpherical) {
    float32_t CosLat, SinLat, rc, xp, zp; /*all local variables */

    /*
     ** Convert geodetic coordinates, (defined by the WGS-84
     ** reference ellipsoid), to Earth Centered Earth Fixed Cartesian
     ** coordinates, and then to spherical coordinates.
     */

    CosLat = arm_cos_f32(DEG2RAD(CoordGeodetic->phi));
    SinLat = arm_sin_f32(DEG2RAD(CoordGeodetic->phi));

    /* compute the local radius of curvature on the WGS-84 reference ellipsoid */
    arm_sqrt_f32(1.0 - Ellip->epssq * SinLat * SinLat, &rc);
    rc = Ellip->a / rc;

    /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */

    xp = (rc + CoordGeodetic->HeightAboveEllipsoid) * CosLat;
    zp = (rc * (1.0 - Ellip->epssq) + CoordGeodetic->HeightAboveEllipsoid) * SinLat;

    /* compute spherical radius and angle lambda and phi of specified point */

    arm_sqrt_f32(xp * xp + zp * zp, &(CoordSpherical->r));
    CoordSpherical->phig = RAD2DEG(asin(zp / CoordSpherical->r)); /* geocentric latitude */
    CoordSpherical->lambda = CoordGeodetic->lambda; /* longitude */
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

static void Geomag_TimelyModifyMagModel(time_t* t, Geomag_MagneticModel_t *MagneticModel, Geomag_MagneticModel_t *TimedMagneticModel) {
    float32_t DecimalYear = Geomag_TimeToYear(t);
    int n, m, index, a, b;
    TimedMagneticModel->EditionDate = MagneticModel->EditionDate;
    TimedMagneticModel->epoch = MagneticModel->epoch;
    TimedMagneticModel->nMax = MagneticModel->nMax;
    TimedMagneticModel->nMaxSecVar = MagneticModel->nMaxSecVar;
    a = TimedMagneticModel->nMaxSecVar;
    b = (a * (a + 1) / 2 + a);
    strcpy(TimedMagneticModel->ModelName, MagneticModel->ModelName);

    for(n = 1; n <= MagneticModel->nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if(index <= b)
            {
                TimedMagneticModel->Main_Field_Coeff_H[index] = MagneticModel->Main_Field_Coeff_H[index] + (DecimalYear - MagneticModel->epoch) * MagneticModel->Secular_Var_Coeff_H[index];
                TimedMagneticModel->Main_Field_Coeff_G[index] = MagneticModel->Main_Field_Coeff_G[index] + (DecimalYear - MagneticModel->epoch) * MagneticModel->Secular_Var_Coeff_G[index];
                TimedMagneticModel->Secular_Var_Coeff_H[index] = MagneticModel->Secular_Var_Coeff_H[index];
                TimedMagneticModel->Secular_Var_Coeff_G[index] = MagneticModel->Secular_Var_Coeff_G[index];
            } else
            {
                TimedMagneticModel->Main_Field_Coeff_H[index] = MagneticModel->Main_Field_Coeff_H[index];
                TimedMagneticModel->Main_Field_Coeff_G[index] = MagneticModel->Main_Field_Coeff_G[index];
            }
        }
    }
}

static void Geomag_ComputeSphericalHarmonicVariables(Geomag_Ellipsoid_t* Ellip, Geomag_CoordSpherical_t* CoordSpherical, int nMax, Geomag_SphericalHarmonicVariables_t* SphVariables) {
    float32_t cos_lambda, sin_lambda;
    int m, n;
    cos_lambda = arm_cos_f32(DEG2RAD(CoordSpherical->lambda));
    sin_lambda = arm_sin_f32(DEG2RAD(CoordSpherical->lambda));
    /* for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2)
    for n  1..nMax-1 (this is much faster than calling pow MAX_N+1 times).      */
    SphVariables->RelativeRadiusPower[0] = (Ellip->re / CoordSpherical->r) * (Ellip->re / CoordSpherical->r);
    for(n = 1; n <= nMax; n++)
    {
        SphVariables->RelativeRadiusPower[n] = SphVariables->RelativeRadiusPower[n - 1] * (Ellip->re / CoordSpherical->r);
    }

    /*
     Compute cos(m*lambda), sin(m*lambda) for m = 0 ... nMax
           cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
           sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
     */
    SphVariables->cos_mlambda[0] = 1.0;
    SphVariables->sin_mlambda[0] = 0.0;

    SphVariables->cos_mlambda[1] = cos_lambda;
    SphVariables->sin_mlambda[1] = sin_lambda;
    for(m = 2; m <= nMax; m++)
    {
        SphVariables->cos_mlambda[m] = SphVariables->cos_mlambda[m - 1] * cos_lambda - SphVariables->sin_mlambda[m - 1] * sin_lambda;
        SphVariables->sin_mlambda[m] = SphVariables->cos_mlambda[m - 1] * sin_lambda + SphVariables->sin_mlambda[m - 1] * cos_lambda;
    }
}

static void Geomag_PcupLow(float32_t *Pcup, float32_t *dPcup, float32_t x, int nMax) {
    int n, m, index, index1, index2, NumTerms;
    float32_t k, z;
    float32_t schmidtQuasiNorm[CALCULATE_NUMTERMS(nMax) + 1];
    Pcup[0] = 1.0;
    dPcup[0] = 0.0;
    /*sin (geocentric latitude) - sin_phi */
    arm_sqrt_f32((1.0 - x) * (1.0 + x), &z);

    /*	 First,	Compute the Gauss-normalized associated Legendre  functions*/
    for(n = 1; n <= nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if(n == m)
            {
                index1 = (n - 1) * n / 2 + m - 1;
                Pcup [index] = z * Pcup[index1];
                dPcup[index] = z * dPcup[index1] + x * Pcup[index1];
            } else if(n == 1 && m == 0)
            {
                index1 = (n - 1) * n / 2 + m;
                Pcup[index] = x * Pcup[index1];
                dPcup[index] = x * dPcup[index1] - z * Pcup[index1];
            } else if(n > 1 && n != m)
            {
                index1 = (n - 2) * (n - 1) / 2 + m;
                index2 = (n - 1) * n / 2 + m;
                if(m > n - 2)
                {
                    Pcup[index] = x * Pcup[index2];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2];
                } else
                {
                    k = (float32_t) (((n - 1) * (n - 1)) - (m * m)) / (float32_t) ((2 * n - 1) * (2 * n - 3));
                    Pcup[index] = x * Pcup[index2] - k * Pcup[index1];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];
                }
            }
        }
    }
    /* Compute the ration between the the Schmidt quasi-normalized associated Legendre
     * functions and the Gauss-normalized version. */

    schmidtQuasiNorm[0] = 1.0;
    for(n = 1; n <= nMax; n++)
    {
        index = (n * (n + 1) / 2);
        index1 = (n - 1) * n / 2;
        /* for m = 0 */
        schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (float32_t) (2 * n - 1) / (float32_t) n;

        for(m = 1; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            index1 = (n * (n + 1) / 2 + m - 1);
            float32_t sqrt;
            arm_sqrt_f32((float32_t) ((n - m + 1) * (m == 1 ? 2 : 1)) / (float32_t) (n + m), &sqrt);
            schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * sqrt;
        }

    }

    /* Converts the  Gauss-normalized associated Legendre
              functions to the Schmidt quasi-normalized version using pre-computed
              relation stored in the variable schmidtQuasiNorm */

    for(n = 1; n <= nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            Pcup[index] = Pcup[index] * schmidtQuasiNorm[index];
            dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index];
            /* The sign is changed since the new WMM routines use derivative with respect to latitude
            insted of co-latitude */
        }
    }
}

static void Geomag_AssociatedLegendreFunction(Geomag_CoordSpherical_t* CoordSpherical, int nMax, Geomag_LegendreFunction_t *LegendreFunction) {
    float32_t sin_phi;
    int FLAG = 1;

    sin_phi = arm_sin_f32(DEG2RAD(CoordSpherical->phig)); /* sin  (geocentric latitude) */

    Geomag_PcupLow(LegendreFunction->Pcup, LegendreFunction->dPcup, sin_phi, nMax);
}

static void Geomag_Summation(Geomag_LegendreFunction_t *LegendreFunction, Geomag_MagneticModel_t *MagneticModel, Geomag_SphericalHarmonicVariables_t* SphVariables, Geomag_CoordSpherical_t* CoordSpherical, Geomag_GeoMagneticElements_t *MagneticResults) {
    int m, n, index;
    float32_t cos_phi;
    MagneticResults->Bz = 0.0;
    MagneticResults->By = 0.0;
    MagneticResults->Bx = 0.0;
    for(n = 1; n <= MagneticModel->nMax; n++)
    {
        for(m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
            /* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
            MagneticResults->Bz -= SphVariables->RelativeRadiusPower[n] *
                    (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->sin_mlambda[m])
                    * (float32_t) (n + 1) * LegendreFunction-> Pcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
            /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
            MagneticResults->By += SphVariables->RelativeRadiusPower[n] *
                    (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->sin_mlambda[m] -
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->cos_mlambda[m])
                    * (float32_t) (m) * LegendreFunction-> Pcup[index];
            /*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
            /* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */

            MagneticResults->Bx -= SphVariables->RelativeRadiusPower[n] *
                    (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->sin_mlambda[m])
                    * LegendreFunction-> dPcup[index];
        }
    }

    cos_phi = arm_cos_f32(DEG2RAD(CoordSpherical->phig));

    MagneticResults->By = MagneticResults->By / cos_phi;

    // TODO: Should we be checking for geographic poles? I do not think it is likely to be exactly +-90 (If it is, it calls one more function)
}

static void Geomag_RotateMagneticVector(Geomag_CoordSpherical_t* CoordSpherical, Geomag_CoordGeodetic_t* CoordGeodetic, Geomag_GeoMagneticElements_t* MagneticResultsSph, Geomag_GeoMagneticElements_t *MagneticResultsGeo) {
    float32_t Psi;
    /* Difference between the spherical and Geodetic latitudes */
    Psi = (PI / 180) * (CoordSpherical->phig - CoordGeodetic->phi);

    /* Rotate spherical field components to the Geodetic system */
    Vec_RotateSpher((Vec3D_t*)&(MagneticResultsSph->Bx), -Psi, 0, (Vec3D_t*)&(MagneticResultsGeo->Bx));
    /* original:
    MagneticResultsGeo->Bz = MagneticResultsSph->Bx * arm_sin_f32(Psi) + MagneticResultsSph->Bz * arm_cos_f32(Psi);
    MagneticResultsGeo->Bx = MagneticResultsSph->Bx * arm_cos_f32(Psi) - MagneticResultsSph->Bz * arm_sin_f32(Psi);
    MagneticResultsGeo->By = MagneticResultsSph->By;*/
}

static void Geomag_CalculateGeomagneticElements(Geomag_GeoMagneticElements_t *MagneticResultsGeo, Geomag_GeoMagneticElements_t *GeoMagneticElements) {
    arm_sqrt_f32(MagneticResultsGeo->Bx * MagneticResultsGeo->Bx + MagneticResultsGeo->By * MagneticResultsGeo->By, &(GeoMagneticElements->H));
    arm_sqrt_f32(GeoMagneticElements->H * GeoMagneticElements->H + MagneticResultsGeo->Bz * MagneticResultsGeo->Bz, &(GeoMagneticElements->F));
    arm_atan2_f32(GeoMagneticElements->By, GeoMagneticElements->Bx, &(GeoMagneticElements->Decl));
    GeoMagneticElements->Decl = RAD2DEG(GeoMagneticElements->Decl);
    arm_atan2_f32(GeoMagneticElements->Bz, GeoMagneticElements->H, &(GeoMagneticElements->Incl));
    GeoMagneticElements->Incl = RAD2DEG(GeoMagneticElements->Incl);
}

float Geomag_SlowFactorial(int n)
{
    int i;

    float res;

    res = 1;
    for (i=2; i<=n; i++)
        res *= i;
    
    return res;
}

float Geomag_CalcLegendere(int n, int m, float x)
{
    float coef;

    coef = 1;
    if(m > 0)
        coef = sqrtf(2 * Geomag_SlowFactorial(n - m) / Geomag_SlowFactorial(n + m));
    

}

float Geomag_ComputePotentialAtPoint(Geomag_Ellipsoid_t* Ellip, Geomag_CoordSpherical_t* CoordSpherical, Geomag_MagneticModel_t* TimedMagneticModel)
{
    int i, n, m;

    float res;

    float altterm, fieldterm, legendreterm;

    res = 0;
    for(i=1, n=1; i<=TimedMagneticModel->nMax; n++)
    {
        for(m=0; m<=n; m++, i++)
        {
            altterm = powf(Ellip->re / CoordSpherical->r, n+1);
            fieldterm = 0;
            fieldterm += TimedMagneticModel->Main_Field_Coeff_G[i] * cosf(CoordSpherical->lambda * m);
            fieldterm += TimedMagneticModel->Main_Field_Coeff_H[i] * sinf(CoordSpherical->lambda * m);
            legendreterm = Geomag_CalcLegendere(n, m, cosf(CoordSpherical->phig));
            res += altterm * fieldterm * legendreterm;
            //printf("i=%d, n=%d, m=%d, g=%f, h=%f\n", i, n, m, TimedMagneticModel->Main_Field_Coeff_G[i], TimedMagneticModel->Main_Field_Coeff_H[i]);
        }
    }

    return res * Ellip->re;
}

static void Geomag_Compute(Geomag_Ellipsoid_t* Ellip, Geomag_CoordSpherical_t* CoordSpherical, Geomag_CoordGeodetic_t* CoordGeodetic, Geomag_MagneticModel_t* TimedMagneticModel, Geomag_GeoMagneticElements_t *GeomagElements)
{   
#if 0
    Geomag_LegendreFunction_t LegendreFunction;
    Geomag_SphericalHarmonicVariables_t SphVariables;
    Geomag_GeoMagneticElements_t MagneticResultsSph, MagneticResultsGeo;

    Geomag_ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, TimedMagneticModel->nMax, &SphVariables); /* Compute Spherical Harmonic variables  */
    Geomag_AssociatedLegendreFunction(CoordSpherical, TimedMagneticModel->nMax, &LegendreFunction); /* Compute ALF  */
    Geomag_Summation(&LegendreFunction, TimedMagneticModel, &SphVariables, CoordSpherical, &MagneticResultsSph); /* Accumulate the spherical harmonic coefficients*/
    Geomag_RotateMagneticVector(CoordSpherical, CoordGeodetic, &MagneticResultsSph, &MagneticResultsGeo); /* Map the computed Magnetic fields to Geodeitic coordinates  */
    Geomag_CalculateGeomagneticElements(&MagneticResultsGeo, GeomagElements); /* Calculate the Geomagnetic elements, Equation 19 , WMM Technical report */
#else
    const float epsilon = 0.01;

    float p, px, py, pz, dx, dy, dz;

    p = Geomag_ComputePotentialAtPoint(Ellip, CoordSpherical, TimedMagneticModel);
    
    CoordSpherical->lambda += epsilon;
    px = Geomag_ComputePotentialAtPoint(Ellip, CoordSpherical, TimedMagneticModel) - p;
    CoordSpherical->lambda -= epsilon;

    CoordSpherical->phig += epsilon;
    py = Geomag_ComputePotentialAtPoint(Ellip, CoordSpherical, TimedMagneticModel) - p;
    CoordSpherical->phig -= epsilon;

    CoordSpherical->r += epsilon;
    pz = Geomag_ComputePotentialAtPoint(Ellip, CoordSpherical, TimedMagneticModel) - p;
    CoordSpherical->r -= epsilon;

    dx = px / epsilon;
    dy = py / epsilon;
    dz = pz / epsilon;

    //printf("p: %f, %f, %f.\n", px, py, pz);
    //printf("d: %f, %f, %f.\n", dx, dy, dz);

    GeomagElements->Bx = dx;
    GeomagElements->By = dy;
    GeomagElements->Bz = dz;

    GeomagElements->F = p;

    Geomag_SphericalToCartesian(&GeomagElements->Bx, &GeomagElements->g);
#endif
}

static void Geomag_HorizontalToEquatorial(MAGtype_GeoMagneticElements* MagHorizontal, Vec3D_t *MagEquatorial, float32_t theta, float32_t phi_p)
{
    float lambda, phi;
    Vec3D_t spherical;

    lambda = DEG2RAD(phi_p);
    phi = DEG2RAD(theta);

    MagEquatorial->X =
        MagHorizontal->X * -sinf(phi) *  cosf(lambda) + 
        MagHorizontal->Y * -sinf(phi) *  sinf(lambda) +
        MagHorizontal->Z *  sinf(phi) *  1.0         ;

    MagEquatorial->Y =
        MagHorizontal->X *  1.0       *  sinf(lambda) + 
        MagHorizontal->Y *  1.0       *  cosf(lambda) +
        MagHorizontal->Z *  0.0       *  0.0         ;

    MagEquatorial->Z =
        MagHorizontal->X * -cosf(phi) *  cosf(lambda) + 
        MagHorizontal->Y * -cosf(phi) *  sinf(lambda) +
        MagHorizontal->Z * -sinf(phi) *  1.0         ; 

    return;

    spherical.X =  MagHorizontal->Y;
    spherical.Y = -MagHorizontal->Z;
    spherical.Z =  MagHorizontal->X;

    Vec_RotateSpher(&spherical, phi, lambda, MagEquatorial);
}

void Geomag_RunTests(const char *filename)
{
    int i;

    FILE *ptr;
    float year, altitude, latitude, longitude, testp;
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
    while (fscanf(ptr, "%f %f %f %f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f\n", 
                  &year, &altitude, &latitude, &longitude, &testp) == 5) 
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
        p = vec.X;

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
