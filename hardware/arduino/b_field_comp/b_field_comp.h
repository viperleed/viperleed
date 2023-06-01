/*
ViPErLEED - Magnetic-field compensation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 16.05.2023
---------------------
*/


#ifndef _VIPERLEED_B_FIELD_COMP
#define _VIPERLEED_B_FIELD_COMP

#include <Arduino.h>
#include <SPI.h>
#include "pwm.h"            // for set_pwm_frequency, set_coil_current
#include "TLE7209.h"        // for motor driver functions

// Include some custom libraries. They are in the ../lib directory.
// This does NOT WORK CORRECTLY when using the Arduino IDE, as the
// latter copies all files into a temporary directory before compiling.
// If you want to use the Arduino IDE, you will have to manually
// Sketch>Add file... and include the following. Bear in mind, however,
// that this DUPLICATES the files in the b_field_comp folder. Thus any
// edit in the original library files WILL NOT BE REFLECTED in the files
// added manually. As an alternative, use upload_sketch.py for compiling,
// and/or uploading, as that one automatically takes care of linking
// library root directories before compilation.
#include "arduino_utils.h"  // for setChipSelectHigh


// Tell the compiler that 'set_signed_pwm_value' is declared in another file
// (i.e., pwm.ino). This is necessary for the CoilClass declaration below
extern byte set_signed_pwm_value(double, byte, byte*);

// Same for the TLE7209 functions used in the driver
extern TLE7209_Error TLE7209readIDandVersion(byte, byte*);
extern TLE7209_Error TLE7209readDiagnosticRegister(byte, byte*);


// The next two pins SHOULD POSSIBLY NOT BE CHANGED. Changing these
// requires picking a different timer/counter module (TC1 or TC3
// instead of the currently used TC4)
#define COIL_1_PWM          6      // PWM output, i.e., voltage value; Also A7;  PD7 on Atmega32U4
#define COIL_1_PWM_REGISTER OCR4D  // WARNING: this is directly related to the choice of COIL_1_PWM

#define COIL_2_PWM          10     // PWM output, i.e., voltage value; Also A10; PB6 on Atmega32U4
#define COIL_2_PWM_REGISTER OCR4B  // WARNING: this is directly related to the choice of COIL_2_PWM

// Current direction: positive or negative?
#define COIL_1_SIGN 18  // Used for INA on shunt; PF7 on ATmega32U4
#define COIL_2_SIGN 19  // Used for INA on shunt; PF6 on ATmega32U4

#define COIL_1_SPI_CS 11     // For SPI communication; PB7 on ATmega32U4
#define COIL_2_SPI_CS 12     // For SPI communication; PD6 on ATmega32U4



class MotorDriver{
    public:
        MotorDriver(byte chip_select_pin) : spi_cs_pin(chip_select_pin) {};
        
        void setup() {
            setChipSelectHigh(spi_cs_pin);
        };
        
        TLE7209_Error get_version(byte* version){
            return TLE7209readIDandVersion(spi_cs_pin, version);
        };
        
        TLE7209_Error get_diagnostic_info(byte* info){
            return TLE7209readDiagnosticRegister(spi_cs_pin, info);
        };

    private:
        const byte spi_cs_pin;
};



class Coil {
    public:
        const MotorDriver driver;

        // Class constructor
        Coil(byte pwm, byte _pwm_register, byte sign, byte spi_cs)
           : driver(spi_cs), pwm_pin(pwm), pwm_register(_pwm_register),
             pwm_sign_pin(sign) {};

        void setup() {
            pinMode(pwm_pin, OUTPUT);
            pinMode(pwm_sign_pin, OUTPUT);
            driver.setup();
        };

        byte set_current(double coil_current) {
            return set_signed_pwm_value(coil_current, pwm_sign_pin,
                                        &pwm_register);
        };

    private:
        // Once the constructor has been called, 'pwm_pin' will be
        // initialized with the desired coil (COIL_1_PWM, COIL_2_PWM)
        const byte pwm_pin;
        const byte pwm_register;
        const byte pwm_sign_pin;
};


// Create two global instances of CoilClass with respective initializers
Coil coil_1(COIL_1_PWM, COIL_1_PWM_REGISTER, COIL_1_SIGN, COIL_1_SPI_CS);
Coil coil_2(COIL_2_PWM, COIL_2_PWM_REGISTER, COIL_2_SIGN, COIL_2_SPI_CS);


#endif   // _VIPERLEED_B_FIELD_COMP
