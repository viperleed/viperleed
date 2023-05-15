/*
Test firmware for Arduino dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/


#include <Arduino.h>  // for interrupts()/noInterrupts()
#include <SPI.h>
#include "pwm.h"
#include "TLE7209.h"

                              // TODO: #include other defines, e.g. 'adc.h', 'struct.h', ...

#define DEBUG true            // Debug mode, writes to serial line, for use in serial monitor


void setup()
{
    #if DEBUG
        Serial.setTimeout(100);
        Serial.begin(9600);           // opens serial port, sets data rate to 9600 bps
    #endif

    set_pwm_frequency(20);            // 20 kHz
    set_coil_current(0.625, COIL_1);
    set_coil_current(0.25,  COIL_2);

    pinMode(COIL_1, OUTPUT);          // Define PD7 (OC4D) as output
    pinMode(COIL_2, OUTPUT);          // Define PB6 (OC4B) as output
    pinMode(COIL_1_SIGN, OUTPUT);     // Define PF7 as output
    pinMode(COIL_2_SIGN, OUTPUT);     // Define PF6 as output

    pinMode(TLE_CHIPSELECT, OUTPUT);
    SPI.begin();                      // Initializes the SPI bus (SCK and MOSI as OUTPUT)
    pinMode(MISO, INPUT);             // MISO = pin PB3
}


void loop()
{
    uint8_t byteRead;
    TLE7209_Error errcode = TLE7209_NoError;
    
    delay(1000);

    #if DEBUG
        Serial.println("Heartbeat\n");
    #endif

    errcode = TLE7209readIDandVersion(TLE_CHIPSELECT, &byteRead);
    delayMicroseconds(50);

    byteRead = TLE7209readDiagnosticRegister(TLE_CHIPSELECT, &byteRead);
    delayMicroseconds(50);
}


byte set_coil_current(double coil_current, uint8_t coil){
    /**Set fraction of maximum current in a coil.

    Parameters
    ----------
    coil_current : double
        Fraction of maximum current. Should be between zero and one.
    coil : {COIL_1, COIL_2}
        Which coil's current should be set.

    Returns
    -------
    error_code : byte
        0 for no error
        2 for coil_current out-of-range
        3 for invalid coil
    **/
    uint8_t *_reg_addr;
    byte sign_select_pin;

    if (coil_current < -1 || coil_current > 1) return 2;                        // TODO: could make these return values into error codes, similar to the driver codes, or use a bunch of defines
    switch(coil)
    {
        case COIL_1:
            _reg_addr = &OCR4D;
            sign_select_pin = COIL_1_SIGN;
            break;
        case COIL_2:
            _reg_addr = &OCR4B;
            sign_select_pin = COIL_2_SIGN;
            break;
        default:
            return 3;
    }

    // Notice that the coil_current, i.e., the time-averaged
    // value of the signal from the PWM, is exactly the same
    // as the duty cycle of the PWM itself.
    // Set PWM duty cycle for channel at `_reg_addr`:
    //    `duty_cycle` == `coil_current` = TC4H:`_register` / TC4H:OCR4C
    set_current_sign(coil_current, sign_select_pin);
    set_ten_bit_value(pwm_clock_divider * coil_current, _reg_addr);
    return 0;
}



