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
// (i.e., pwm.ino). This is necessary for the 'Coil' class declaration below
extern byte set_signed_pwm_value(double, byte, byte*);
extern void enable_pwm_channel(TC4_PWM_CHANNEL, bool);

// Same for the TLE7209 functions used in the driver
extern void TLE7209reset(byte);
extern TLE7209_Error TLE7209readIDandVersion(byte, byte *);
extern TLE7209_Error TLE7209readDiagnosticRegister(byte, byte *);


// The pins belonging to 'COIL_1_PWM' and 'COIL_2_PWM' should not be changed.
// Changing these may require choosing a different Timer/Counter module, e.g.
// TC1 or TC3 instead of the currently used TC4.
#define COIL_1_PWM               6   // PWM output 1, i.e., voltage value; Also A7;  PD7 on Atmega32U4                          //  DO NOT CHANGE (OC4D PWM output)
#define COIL_1_DISABLE          21   // Could later on be an alias of signal "COIL_1_SIGN"
#define COIL_1_PWM_REGISTER     FAST_PWM_CH_1_REG  // WARNING: this is directly related to the choice of COIL_1_PWM

#define COIL_2_PWM              10   // PWM output 1, i.e., voltage value; Also A10; PB6 on Atmega32U4                          //  DO NOT CHANGE (OC4B PWM output)
#define COIL_2_DISABLE          22   // Could later on be an alias of signal "COIL_2_SIGN"
#define COIL_2_PWM_REGISTER     FAST_PWM_CH_2_REG  // WARNING: this is directly related to the choice of COIL_2_PWM

// Current direction: positive or negative?
#define COIL_1_SIGN             18  // Used for INA on shunt; PF7 on ATmega32U4
#define COIL_2_SIGN             19  // Used for INA on shunt; PF6 on ATmega32U4

#define COIL_1_SPI_CS           11  // For SPI communication; PB7 on ATmega32U4
#define COIL_2_SPI_CS           12  // For SPI communication; PD6 on ATmega32U4



class MotorDriver{
    public:
        MotorDriver(byte chip_select_pin, byte _disable_pin)
        : spi_cs_pin(chip_select_pin), disable_pin(_disable_pin) {};

        void setup() {
            setChipSelectHigh(spi_cs_pin);
            pinMode(disable_pin, OUTPUT);    // Use for TLE7209 disable pin 'DIS'
            digitalWrite(disable_pin, LOW);  // 'DIS' == 0: Activate TLE7209 output drivers
        };

        void reset() {
            TLE7209reset(disable_pin);
        };

        TLE7209_Error get_version(byte* version) {
            return TLE7209readIDandVersion(spi_cs_pin, version);
        };

        TLE7209_Error get_diagnostic_info(byte* info) {
            return TLE7209readDiagnosticRegister(spi_cs_pin, info);
        };

    private:
        const byte spi_cs_pin;
        const byte disable_pin;
};


class Coil {
    public:
        const MotorDriver driver;

        // Class constructor, including member initializer list
        Coil(byte pwm, byte sign, byte spi_cs, byte disable,
             byte _pwm_register_addr)
           : driver(spi_cs, disable), pwm_pin(pwm),
             pwm_register_addr(_pwm_register_addr),
             pwm_sign_pin(sign) {};

        void setup() {
            pinMode(pwm_pin, OUTPUT);
            pinMode(pwm_sign_pin, OUTPUT);

            // Reset the TLE7209 to clear any previous error condition
            driver.reset();
            driver.setup();
            set_current(0.0);            
            enable_pwm_channel(tc4_pwm_channel, true);      
        };

        byte set_current(double coil_current) {
            byte err = set_signed_pwm_value(coil_current, pwm_sign_pin,
                                        tc4_reg_addr);
            if(err)
                return err;
            last_current_setpoint = coil_current;
        };

        double get_current() {
            return last_current_setpoint; 
        }
    private:
        // Once the constructor has been called, 'pwm_pin' will be
        // initialized with the desired coil (COIL_1_PWM, COIL_2_PWM)
        const byte pwm_pin;
        const byte pwm_register_addr;
        const byte pwm_sign_pin;
        double last_current_setpoint;
};


// Create two global instances of class 'Coil' with respective initializers
Coil coil_1(COIL_1_PWM, COIL_1_SIGN, COIL_1_SPI_CS, COIL_1_DISABLE, COIL_1_PWM_REGISTER);
Coil coil_2(COIL_2_PWM, COIL_2_SIGN, COIL_2_SPI_CS, COIL_2_DISABLE, COIL_2_PWM_REGISTER);


#endif   // _VIPERLEED_B_FIELD_COMP
