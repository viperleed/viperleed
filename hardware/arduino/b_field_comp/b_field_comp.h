/*
ViPErLEED - Magnetic-field compensation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 16.05.2023
---------------------
*/


#ifndef _VIPERLEED_B_FIELD_COMP
#define _VIPERLEED_B_FIELD_COMP

#include <Arduino.h>       // for interrupts()/noInterrupts()
#include <SPI.h>
#include "pwm.h"           // for set_pwm_frequency, set_coil_current
#include "TLE7209.h"

// TODO: move here all the coil-related defines, and make
// them into a class that also contains relevant methods
#define TLE_CHIPSELECT 11     // For SPI communication; Also PB7                // TODO: need two of these; both inside the coil class; name suggestion: driver_spi or something along these lines?



#endif
