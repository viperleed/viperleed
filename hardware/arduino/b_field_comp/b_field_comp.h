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


// The next two pins SHOULD POSSIBLY NOT BE CHANGED. Changing these
// requires picking a different timer/counter module (currently TC4
// TC1 or TC3)
#define COIL_1  6   // PWM output, i.e., voltage value; Also A7;  PD7 on Atmega32U4
#define COIL_2 10   // PWM output, i.e., voltage value; Also A10; PB6 on Atmega32U4

// Current direction: positive or negative?
#define COIL_1_SIGN 18  // Used for INA on shunt; PF7 on ATmega32U4 
#define COIL_2_SIGN 19  // Used for INA on shunt; PF6 on ATmega32U4

#define TLE_CHIPSELECT_1 11     // For SPI communication; PB7 on ATmega32U4 
#define TLE_CHIPSELECT_2 12     // For SPI communication; PD6 on ATmega32U4 
#endif
