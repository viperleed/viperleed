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



class CoilClass {
public:
  // Once the constructor has been called, 'coil' will be initialized with
  // the desired coil (COIL_1, COIL_2)
  byte coil;
  
  // Class constructor
  CoilClass(byte _coil) {
    coil = _coil;
  }

  void setup() {
    pinMode(coil, OUTPUT);
  };

  byte set_current(double coil_current) {
    return set_coil_current(coil_current, coil);
  };
};

// Create two global instances of CoilClass with respective initializers
CoilClass Coil1(COIL_1), Coil2(COIL_2);


#endif
