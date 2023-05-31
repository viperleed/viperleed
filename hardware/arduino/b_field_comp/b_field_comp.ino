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

    Coil1.set_current(0.625);
    Coil2.set_current(0.25);

    Coil1.setup();                    // Define PD7 and PB6 as output (used for PWM)
    Coil2.setup();

    pinMode(COIL_1_SIGN, OUTPUT);     // Define PF7 as output
    pinMode(COIL_2_SIGN, OUTPUT);     // Define PF6 as output

    //setChipSelectHigh(TLE_CHIPSELECT);                                          // TODO: The compiler cannot find 'arduino_utils.h'
    digitalWrite(TLE_CHIPSELECT_1, HIGH);
    digitalWrite(TLE_CHIPSELECT_2, HIGH);

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

    errcode = TLE7209readIDandVersion(TLE_CHIPSELECT_1, &byteRead);
    delayMicroseconds(50);

    byteRead = TLE7209readDiagnosticRegister(TLE_CHIPSELECT_1, &byteRead);
    delayMicroseconds(50);


    /*
    TODO:
        - Call 'set_coil_current' and adjust current as needed
        - Check return code, transmit error message over serial console
        - Integrate ADC readout from INA214 and current shunt
    */
    // set_coil_current(double coil_current, uint8_t coil)
}




