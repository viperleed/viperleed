/*
Test firmware for Arduino dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/


#include <Arduino.h>  // for interrupts()/noInterrupts()
#include <SPI.h>
#include "pwm.h"      // for set_pwm_frequency, set_coil_current
#include "TLE7209.h"

// TODO: #include other defines, e.g. 'adc.h', 'struct.h', ...

#define DEBUG true    // Debug mode, writes to serial line, for use in serial monitor


void setup()
{
    #if DEBUG
        Serial.setTimeout(100);
        Serial.begin(9600);           // opens serial port, sets data rate to 9600 bps
    #endif

    set_pwm_frequency(20);  // kHz
    set_coil_current(0.625, COIL_1);                                            // TODO: could become a .set_current method of coil_t objects? Similarly, the two pinMode calls below could be grouped into a .setup method the coil_t object
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


    /*
    TODO:
        - Call 'set_coil_current' and adjust current as needed
        - Check return code, transmit error message over serial console
        - Integrate ADC readout from INA214 and current shunt
    */
    // set_coil_current(double coil_current, uint8_t coil)
}




