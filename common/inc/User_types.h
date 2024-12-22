#ifndef _user_types_h
#define _user_types_h

#ifdef __cplusplus
extern "C"
{
#endif

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* INCLUDES
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
#include "arm_math.h"
#include <stdint.h>

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* EXTERNAL DEFINES
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/* NULL_PTR define with a void pointer to zero (the usage does not lead to compile warnings about implicit casts) */
# ifndef NULL_PTR
#  define NULL_PTR ((void *)0)
# endif

#ifndef PI
#  define PI ((float)3.141592653589793)
#endif

/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* EXTERNAL TYPES DECLARATIONS
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
typedef enum MagnetometersType {
  Magnitometer_0,
  Magnitometer_1,
  Magnitometer_2,
  NO_OF_MMETER_IDS
} MagnetometersType;

typedef enum ES_ReturnType {
  E_OK,
  E_PENDING,
  E_UNKNOWN_FAILURE,
  E_NOT_INIT,
  E_INV_ARG,
  E_INV_CFG,
  E_ERROR,
  E_MMETER_READING,
  E_FILTERING,
  E_DATA_NORMALIZE,
  E_BDOT_CALCULATING,
  E_MTORQ_UPDATE,
  E_NOT_OK
} ES_ReturnType;

typedef struct Magnetorquer_Axis_t {
  double AXIS_X;
  double AXIS_Y;
  double AXIS_Z;  
} Magnetorquer_Axis_t;

typedef struct Compass_Axis_t{
  double AXIS_X;
  double AXIS_Y;
  double AXIS_Z;
} Compass_Axis_t;

typedef enum Compass_mode_t{
  E_CONTINIOUS,
  E_SINGLE,
  E_IDLE,
  E_SLEEP
} Compass_mode_t;

typedef enum Compass_Measuring_mode_t{
  E_MASK,
  E_NORMAL,
  E_POSITIVE_BIAS,
  E_NEGATIVE_BIAS
} Compass_Measuring_mode_t;

typedef enum Compass_ID_Def_t{
  E_A_REG,
  E_B_REG,
  E_C_REG,
  E_COMBINED
} Compass_ID_Def_t;

typedef enum Compass_Set_Range_t{
	E_90uT,
	E_130uT,
	E_190uT,
	E_250uT,
	E_400uT,
	E_470uT,
	E_560uT,
	E_810uT
}Compass_Set_Range_t;

typedef enum Compass_Set_Bandwidth_t {
	E_0_75HZ,  /* 0.75 Hz  XXX called "1" for now XXX */
	E_1_5HZ,   /* 1.5 Hz  XXX called "2" for now XXX */
	E_3HZ,
	E_7_5HZ,   /* 7.5 Hz  XXX called "8" for now XXX */
	E_15HZ,
	E_30HZ,
	E_75HZ  
}Compass_Set_Bandwidth_t;

typedef enum Compass_Samples_Averaged_t
{
    E_MEAS_AVG,
    E_MEAS_AVG_1,
    E_MEAS_AVG_2,
    E_MEAS_AVG_4,
    E_MEAS_AVG_8
} Compass_Samples_Averaged_t;        

typedef enum MagnetorquerType
{
    Magnetorquer_0,
    Magnetorquer_1,
    Magnetorquer_2,
    NO_OF_MTORQ_IDS
} MagnetorquerType;

typedef enum
{
  SEN_SUCCESS =0x01,
  SEN_ERROR   =0x00,
} status_t;

typedef struct
{
  int16_t AXIS_X;
  int16_t AXIS_Y;
  int16_t AXIS_Z;
} AxesRaw_t;

typedef struct
{
  int16_t Temp_X;
  int16_t Temp_Y;
  int16_t Temp_Z;
} Temperature_t;

typedef union Vec3D_t
{
    struct
    {
        float32_t X, Y, Z;
	};
    float32_t Vec[3];
} Vec3D_t;

#ifdef __cplusplus
}
#endif

#endif
