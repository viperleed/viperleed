/*
Test firmware for Arduino Micro dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/


#include "b_field_comp.h"             // Globals, #defines, etc.
//#include "arduino_utils.h"          // from '../lib'; for setChipSelectHigh

#define DEBUG true                    // Debug mode, writes to serial line, for use in serial monitor


void setup()
{
    #if DEBUG
        Serial.setTimeout(100);
        Serial.begin(9600);           // Open serial port, set data rate to 9600 bps
        delay(2000);                  // The IDE needs some time to connect to the Serial Monitor
    #endif

    set_pwm_frequency(20000);
    enable_fast_pwm_mode();
    enable_pwm_channel(TC4_PWM_CH_D, true);         // Enable PWM on OC4D/Pin 6 ('COIL_1_PWM')
    enable_pwm_channel(TC4_PWM_CH_B, true);         // Enable PWM on OC4B/Pin 10 ('COIL_2_PWM')    

    coil_1.setup();                                 // Define PD7 and PB6 as output (used for PWM)
    coil_2.setup();                                 // The setup method also sets the current to 0

    coil_1.set_current(0.625);                      // Make sure to call .set_current *after* .setup
    coil_2.set_current(0.25);

    SPI.begin();                                    // Initialize the SPI bus (SCK and MOSI as OUTPUT)
    pinMode(MISO, INPUT);                           // MISO = pin PB3
}


void loop()
{
    uint8_t byteRead;
    TLE7209_Error errcode = TLE7209_NoError;

    delay(1000);

    #if DEBUG
        Serial.println("Heartbeat\n");
    #endif

    errcode = coil_2.driver.get_version(&byteRead);
    delayMicroseconds(50);

    byteRead = coil_2.driver.get_diagnostic_info(&byteRead);
    delayMicroseconds(50);


    /*
    TODO:
        - Call '.set_current' method and adjust current as needed
        - Check return code, transmit error message over serial console
        - Integrate ADC readout from INA214 and current shunt
        - Periodically monitor the DIA_REG for short-circuit, overcurrent/-temperature
          If a fault occurs, call driver.reset() to clear the error condition and re-enable the output drivers

    */
    
}




