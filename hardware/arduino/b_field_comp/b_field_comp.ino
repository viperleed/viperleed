/*
Test firmware for Arduino Micro dual channel PWM generation
---------------------
Authors: Michele Riva, Christoph Pfungen, Stefan Mitterhöfer, Florian Dörr, Tun Sinner
15.05.2023, MR, CP: First version
04.12.2024, SM: Modify setup to changes in class coil, include measurement of characteristic curve
01.04.2025 FD, TS: Addition of a state machine + communication protocol
---------------------
*/

#include <SPI.h>

#define DEBUG false                  // Debug mode, writes to serial line, for use in serial monitor

#include "states-def.h"				 // Basic state machine definitions
#include "viper-serial.h"			 // Serial communication functions
#include "arduino_utils.h"           // from '../lib'; for setChipSelectHigh
#include "b_field_comp.h"            // Globals, #defines, etc.


void setup()
{
  Serial.setTimeout(100);
  Serial.begin(9600);          // Open serial port, set data rate to 9600 bps
  delay(2000);                 // The IDE needs some time to connect to the Serial Monitor

  /*

  //analogReference(INTERNAL);       // Enable 2.56V bandgap reference on Aref pin

  */

  coil_1.setup();
  coil_2.setup();

  // coil_1.set_duty_cycle(-0.1);
  // coil_2.set_duty_cycle(0.2);

  delay(100);

}

void loop()
{
    //uint8_t byteRead;
    // TLE7209_Error errcode = TLE7209_NoError;
    float current, temperature, shuntvoltage, voltage;							// TODO: globals?
    float setpoint;

    // delay(1000);

    #if DEBUG
        Serial.println("Heartbeat\n");
    #endif


	/***
    setpoint = -1.000;
    Serial.println("duty cycle;voltage;current");
    do {
      coil_1.set_duty_cycle(setpoint);
      delay(100);
      voltage = coil_1.get_voltage();
      current = coil_1.get_current();

      Serial.print(setpoint,3);
      Serial.print(";");
      Serial.print(voltage);
      Serial.print(";");
      Serial.println(current,6);

      setpoint += 0.001;
    }while(setpoint < 1.0);
	***/




	readFromSerial();
	updateState();
	compensateForSupplyVoltage();													// TODO: Really here?
	compensateForThermalDrift();





    // do {
    //   coil_1.set_duty_cycle(setpoint);
    //   delay(100);
    //   voltage = coil_1.get_voltage();
    //   current = coil_1.get_current();

    //   Serial.print(setpoint,3);
    //   Serial.print(";");
    //   Serial.print(voltage);
    //   Serial.print(";");
    //   Serial.println(current,6);

    //   setpoint -= 0.01;
    // }while(setpoint >= -1);

    // while(1) {
    //   setpoint = -0.2;
    //   coil_1.set_duty_cycle(setpoint);
    //   coil_2.set_duty_cycle(setpoint);
    //   Serial.print("setpoint = ");
    //   Serial.println(setpoint);
    //   delay(100);
    //   voltage = coil_1.get_voltage();
    //   current = coil_1.get_current();
    //   Serial.print("voltage = ");
    //   Serial.print(voltage);
    //   Serial.println(" V");
    //   Serial.print("current = ");
    //   Serial.print(current, 6);
    //   Serial.println(" A");
    //   setpoint = setpoint * 2;
    //   coil_1.set_duty_cycle(-setpoint);
    //   coil_2.set_duty_cycle(-setpoint);
    //   Serial.print("setpoint = ");
    //   Serial.println(-setpoint);
    //   delay(100);
    //   voltage = coil_1.get_voltage();
    //   current = coil_1.get_current();
    //   Serial.print("voltage = ");
    //   Serial.print(voltage);
    //   Serial.println(" V");
    //   Serial.print("current = ");
    //   Serial.print(current, 6);
    //   Serial.println(" A");
    // }





	/***
    // shuntvoltage = coil_1.measurement.get_shunt_voltage();
    // current = coil_1.get_current();
    // temperature = coil_1.get_ambient_temperature();
    // voltage = coil_1.get_voltage();

    // #if DEBUG
        // Serial.print("shuntvoltage = ");
        // Serial.print(shuntvoltage);
        // Serial.println(" mV");
        // Serial.print("current = ");
        // Serial.print(current, 6);
        // Serial.println(" A");
        // Serial.print("temperature = ");
        // Serial.print(temperature);
        // Serial.println(" °C");
        // Serial.print("voltage = ");
        // Serial.print(voltage);
        // Serial.println(" V");
    // #endif
	***/





    /*
    TODO:
        - Call '.set_current' method and adjust current as needed
        - Check return code, transmit error message over serial console
        - Periodically monitor the DIA_REG for short-circuit, overcurrent/-temperature
          If a fault occurs, call driver.reset() to clear the error condition and re-enable the output drivers

    */

}



bool isAllowedCommand() {
	return true;																// TODO: check that the received command is one of the commands that initiates one of the defined states
}

void updateState() {
	if (not newMessage) return; // No new message, therefore do not set state

	if (msgLength > 1) return; // We received data, not a command

	switch(data_received[0]){
		case PC_SOMETHING:													// TODO: fill in commands for checks and set states
			doSomething();
			break;
		case PC_SOMETHING_ELSE:
			somethingElse();
			currentState = SOME_STATE;
			break;
		case
	}

	newMessage = false;
}

void compensateForSupplyVoltage(){
	/** Measure supply voltage and adjust duty cycle according to it. **/
	// If target current:
	// output_voltage = duty_cycle * supply_voltage
}

void compensateForThermalDrift(){
	/** Compensate for the increased resistance when coils heat up. **/
	// If target current:
	// Only if enough time has passed/not too much current change happened
	// Is change small enough -> calculate resistance
}
