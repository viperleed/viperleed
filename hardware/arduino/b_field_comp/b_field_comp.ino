/*
Test firmware for Arduino dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/


#include "b_field_comp.h"  // globals, #defines, etc.
//#include "arduino_utils.h" // from '../lib'; for setChipSelectHigh

// TODO: #include other defines, e.g. 'adc.h', 'struct.h', ...

#define DEBUG true    // Debug mode, writes to serial line, for use in serial monitor



void setup()
{
    #if DEBUG
        Serial.setTimeout(100);
        Serial.begin(9600);           // opens serial port, sets data rate to 9600 bps
    #endif

    set_pwm_frequency(20);            // kHz

    coil_1.set_current(0.625);
    coil_2.set_current(0.25);

    coil_1.setup();                    // Define PD7 and PB6 as output (used for PWM)
    coil_2.setup();

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

    errcode = coil_1.driver.get_version(&byteRead);
    delayMicroseconds(50);

    byteRead = coil_1.driver.get_diagnostic_info(&byteRead);
    delayMicroseconds(50);


    /*
    TODO:
        - Call '.set_current' method and adjust current as needed
        - Check return code, transmit error message over serial console
        - Integrate ADC readout from INA214 and current shunt
    */
    // set_coil_current(double coil_current, uint8_t coil)
}




