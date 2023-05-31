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


// Tell the compiler that 'set_coil_current' is declared in another file
// This is necessary for the CoilClass declaration below
extern byte set_coil_current(double, uint8_t);


#endif
