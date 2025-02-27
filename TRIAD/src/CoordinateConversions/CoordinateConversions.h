#ifndef _coordinate_conversions_h
#define _coordinate_conversions_h

#include <time.h>

#include <User_types.h>

/* UTC time at which 0deg N and 0deg E was along the X axis in equatorial coordinates.
   This happens at the time of a vernal (spring) equinox.
   (Actually this also happens every day whenever RA 0 crosses the meridian at 0deg N 0deg E,
   but an equinox is a more convenient reference point.) */
#define EQUINOX_TIME ((time_t) 1695451800)

#define EARTH_ROT_SPEED ((float)7.2921159E-5)

#define RAD2DEG(rad)    ((rad)*(180.0L/PI))
#define DEG2RAD(deg)    ((deg)*(PI/180.0L))

void Coord_NEDToECEF(Vec3D_t* NED, Vec3D_t *ECEF, float32_t lat, float32_t lon);
void Coord_ECEFToECI(time_t t, const Vec3D_t* ECEF, Vec3D_t* ECI);
void Coord_ECIToECEF(time_t t, const Vec3D_t* ECI, Vec3D_t* ECEF);

void Coord_TestConversions(void);

#endif