/*
Test firmware for Arduino Micro dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen, Stefan Mitterhöfer
Date: 15.05.2023
---------------------
*/


#include "b_field_comp.h"            // Globals, #defines, etc.
#include "arduino_utils.h"           // from '../lib'; for setChipSelectHigh
//#include "pwm.h"                     // For tc4_sfr_reset, set_pwm_frequency

#define DEBUG true                   // Debug mode, writes to serial line, for use in serial monitor


void setup()
{
    #if DEBUG
        Serial.setTimeout(100);
        Serial.begin(9600);          // Open serial port, set data rate to 9600 bps
        delay(2000);                 // The IDE needs some time to connect to the Serial Monitor
    #endif
    
    setChipSelectHigh(INA_1_SPI_CS); // Initialize CS PIN
    SPI.begin();                     // Initialize the SPI bus (SCK and MOSI as OUTPUT)
    pinMode(MISO, INPUT);            // MISO = pin PB3

    INA_1.reset();
    delayMicroseconds(50);
    INA_1.setup();

    tc4_sfr_reset();                 // Make sure the TC4 registers are reset on startup
    set_pwm_frequency(2000);        // Also enable Fast PWM mode on Timer/Counter4

    coil_1.setup();                  // Define PWM output pins, enable PWM channels
    //coil_2.setup();                  // The setup method also sets the current to 0

    coil_1.set_current(0.0);       // Make sure to call .set_current *after* .setup
    //coil_2.set_current(0.25);


    //analogReference(INTERNAL);       // Enable 2.56V bandgap reference on Aref pin
}


void loop()
{
    uint8_t byteRead;
    uint32_t adc_median_val = 0;
    TLE7209_Error errcode = TLE7209_NoError;
    float current, temperature, shuntvoltage, busvoltage;
    float setpoint;

    delay(1000);

    #if DEBUG
        Serial.println("Heartbeat\n");
    #endif

//    errcode = coil_2.driver.get_version(&byteRead);
//    delayMicroseconds(50);

//    byteRead = coil_2.driver.get_diagnostic_info(&byteRead);
//    delayMicroseconds(50);
    
    setpoint = -1.000;
    Serial.println("duty cycle;shuntvoltage;current");
    do {
      coil_1.set_current(setpoint);
      delay(100);
      shuntvoltage = INA_1.get_voltage_shunt();
      current = INA_1.get_current();
      
      Serial.print(setpoint,3);
      Serial.print(";");
      Serial.print(shuntvoltage);
      Serial.print(";");
      Serial.println(current,6);

      setpoint += 0.001;
    } while(setpoint < 1);

    do {
      coil_1.set_current(setpoint);
      delay(100);
      shuntvoltage = INA_1.get_voltage_shunt();
      current = INA_1.get_current();
      
      Serial.print(setpoint,3);
      Serial.print(";");
      Serial.print(shuntvoltage);
      Serial.print(";");
      Serial.println(current,6);

      setpoint -= 0.001;
    } while(setpoint > -1);
    
    shuntvoltage = INA_1.get_voltage_shunt();
    current = INA_1.get_current();
    temperature = INA_1.get_temp();
    busvoltage = INA_1.get_voltage_bus();

    #if DEBUG
        Serial.print("shuntvoltage = ");
        Serial.print(shuntvoltage);
        Serial.println(" mV");
        Serial.print("current = ");
        Serial.print(current, 6);
        Serial.println(" A");
        Serial.print("temperature = ");
        Serial.print(temperature);
        Serial.println(" °C");
        Serial.print("busvoltage = ");
        Serial.print(busvoltage);
        Serial.println(" V");
    #endif


    /*
    TODO:
        - Call '.set_current' method and adjust current as needed
        - Check return code, transmit error message over serial console
        - Integrate ADC readout from INA214 and current shunt
        - Periodically monitor the DIA_REG for short-circuit, overcurrent/-temperature
          If a fault occurs, call driver.reset() to clear the error condition and re-enable the output drivers

    */
    
}
