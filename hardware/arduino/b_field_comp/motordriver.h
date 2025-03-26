/*
ViPErLEED - Driver for motor driver.
---------------------
Authors: Stefan Mitterhöfer
26.11.2024, SM: First version
---------------------
*/

#ifndef _VIPERLEED_MOTORDRIVER
#define _VIPERLEED_MOTORDRIVER

// Decide here, which Motordriver is chosen: TLE7209, DRV8842
// For each motordriver a subclass is implemented (in separate headerfile)
#define USE_MD_TLE7209    true       // Otherwise DRV8842 is used

#define MD_PWM_PULSE_MIN  10e-6      // 10 µs; measured minimum PWM pulse width
                                     // which can be provided by the motordriver
#define MD_PWM_FREQUENCY_MAX (1/MD_PWM_PULSE_MIN)

#define F_PWM             2e4        // 20 kHz; usual PWM frequency                            
// For duty cycle values, which will generate shorter PWM pulses than 
// 'MD_PWM_PULSE_MIN', the PWM frequency will be decreased automatically to
// avoid to high discrepancy in the set and actual value of the current. 
// This has to be done because of limitations of the motordriver.
#define DUTY_CYCLE_MIN (MD_PWM_PULSE_MIN * F_PWM)    // 2e4 * 10e-6 = 0.2
#define DUTY_CYCLE_MAX (1 - DUTY_CYCLE_MIN)


class MotorDriver {
    public:

        MotorDriver(byte enable_pin)
        : _ENABLE_PIN(enable_pin) {}
        
        MotorDriver(byte chip_select_pin, byte enable_pin)
        : _SPI_CS_PIN(chip_select_pin), _ENABLE_PIN(enable_pin) {}

        virtual void setup() {
            pinMode(_ENABLE_PIN, OUTPUT);    // Use for TLE7209 enable pin 'EN'
            digitalWrite(_ENABLE_PIN, HIGH); // 'DIS' == 0: Activate TLE7209
                                             // output drivers
        }

        void reset() {
            reset(_ENABLE_PIN);
        }

    protected:
        const byte _SPI_CS_PIN = 0;
        const byte _ENABLE_PIN;
        
        /** Reset the TLE7209 after a fault condition occurred.**/
        void reset(byte enable_pin) {
            digitalWrite(_ENABLE_PIN, HIGH);
            delayMicroseconds(20);
            digitalWrite(_ENABLE_PIN, LOW);
            delayMicroseconds(20);
            digitalWrite(_ENABLE_PIN, HIGH);
        }

};

#if USE_MD_TLE7209
  #define MOTORDRIVER  TLE7209
  #include "TLE7209.h"    // for motor driver functions
#else
  #define MOTORDRIVER  DRV8842
  #include "DRV8842.h"    // for motor driver functions
//  #define MOTORDRIVER  MotorDriver
#endif
                                                                                

#endif  // _VIPERLEED_MOTORDRIVER