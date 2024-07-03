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
#include "INA239.h"                  // For configuration and functions of current measurement device INA239

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


// Generic error codes
enum error_t {
    NoError,
    OutOfRange,
    InvalidIOPin,
    InvalidChannel,
    InvalidPrescaler,
    NotImplemented,
};

// Forward declarations
error_t set_signed_pwm_value(double, byte, byte *);
error_t enable_pwm_channel(TC4_PWM_CHANNEL, bool);

void TLE7209reset(byte);
TLE7209_Error TLE7209readIDandVersion(byte, byte *);
TLE7209_Error TLE7209readDiagnosticRegister(byte, byte *);

uint8_t pin_to_tc4_reg_addr(uint8_t);
TC4_PWM_CHANNEL pin_to_tc4_channel(uint8_t);


// The pins belonging to 'COIL_1_PWM' and 'COIL_2_PWM' should not be changed.
// Changing these may require choosing a different Timer/Counter module, e.g.
// TC1 or TC3 instead of the currently used TC4.
#define COIL_1_PWM              10    // PWM output 1, i.e., voltage value; Also A10; PB6 on Atmega32U4                          //  DO NOT CHANGE (OC1B PWM output)
#define COIL_1_DISABLE          2    // Could later on be an alias of signal "COIL_1_SIGN"

#define COIL_2_PWM              10   // PWM output 1, i.e., voltage value; Also A10; PB6 on Atmega32U4                          //  DO NOT CHANGE (OC4B PWM output)
#define COIL_2_DISABLE          3    // Could later on be an alias of signal "COIL_2_SIGN"

// Current direction: positive or negative?
#define COIL_1_SIGN             4    // Used for motor driver; PD4 on ATmega32U4
#define COIL_2_SIGN             5    // Used for motor driver; PC6 on ATmega32U4

#define COIL_1_SPI_CS           11   // For SPI communication; PB7 on ATmega32U4
#define COIL_2_SPI_CS           12   // For SPI communication; PD6 on ATmega32U4


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

                                                                                // TODO: Decide here, which Motordriver is chosen, i.e. TLE209, DRV8842
                                                                                // TODO: Implement in Motordriver Subclasses for each Motordriver
class Coil {
    public:
        const MotorDriver driver;

        // Class constructor, including member initializer list
        // Assign 'tc4_reg_addr', 'tc4_pwm_channel' at the earliest possible time
        Coil(byte pwm, byte sign, byte spi_cs, byte disable)             
           : driver(spi_cs, disable), pwm_pin(pwm),
             pwm_sign_pin(sign) {
            tc4_reg_addr = pin_to_tc4_reg_addr(pwm_pin);
            tc4_pwm_channel = pin_to_tc4_channel(pwm_pin);
        };

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

        double get_current() {                                                  // TODO: Return measured value from INA (no avg, just median)
            
        }
    private:
        byte tc4_reg_addr;
        TC4_PWM_CHANNEL tc4_pwm_channel;       
        const byte pwm_pin;
        const byte pwm_sign_pin;
        double last_current_setpoint;
};


// Create two global instances of class 'Coil' with respective initializers
Coil coil_1(COIL_1_PWM, COIL_1_SIGN, COIL_1_SPI_CS, COIL_1_DISABLE);
Coil coil_2(COIL_2_PWM, COIL_2_SIGN, COIL_2_SPI_CS, COIL_2_DISABLE);


#endif   // _VIPERLEED_B_FIELD_COMP
