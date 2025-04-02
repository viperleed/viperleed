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

// The box ID is an indentifier that is necessary for the PC to know what
// type of Arduino it is handling. 2 is the identifier of a
// b-field-compensation controller. The box ID must never be changed!
#define BOX_ID  2

// Firmware version (MAX: v255.255). CURENTLY: v0.1
#define FIRMWARE_VERSION_MAJOR    0  // max 255
#define FIRMWARE_VERSION_MINOR    1  // max 255



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



    /**Read from serial line and handle the state machine.**/
    // See if there is anything coming from the PC.
	readFromSerial();

	// Then see if the current state needs updating,
	// depending on the message receive
	updateState();

	compensateForSupplyVoltage();													// TODO: Really here?
	compensateForThermalDrift();

    // Finally do what's needed in the currentState
    switch (currentState){
        case STATE_GET_CONFIGURATION:
            getConfiguration();
            break;
		case STATE_SET_SERIAL_NR:
            setSerialNr();
            break;
	}



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
	/**
	Check if the received command is among the commands the arduino can
	process.

	Returns
    -------
    True if the message is acceptable
	**/
    // Check that it is one of the understandable commands
    switch(data_received[0]){
        case PC_CONFIGURATION: break;
        case PC_SET_SERIAL_NR: break;
        default:
            raise(ERROR_MSG_UNKNOWN);
            return false;
	}
    return true;																// TODO: check that the received command is one of the commands that initiates one of the defined states
}


void updateState() {
	/**
    If there is an unprocessed serial message, decide whether
    this requires us to change the state of the Arduino, and
    do some preparation that may be needed before entering
    the new state. All data messages have to be at least two
    bytes long. Messages that are one byte long are assumed
    to be commands.

    Reads
    -----
    newMessage, data_received

    Writes
    ------
    initialTime, currentState, waitingForDataFromPC, calibrationGain

    Goes to state
    -------------
    May go anywhere
    **/
	if (not newMessage) return; // No new message, therefore do not set state

	if (msgLength > 1) return; // We received data, not a command

	switch(data_received[0]){
		case PC_CONFIGURATION:
			initialTime = millis();
			currentState = STATE_GET_CONFIGURATION;
			break;
		case PC_SET_SERIAL_NR:
			waitingForDataFromPC = true;
            initialTime = millis();
			currentState = STATE_SET_SERIAL_NR;
			break;
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


/** Handler of STATE_GET_CONFIGURATION */
void getConfiguration(){
    /**Send box ID, firmware version and serial nr. to PC.
    The serial number is read from the EEPROM.

    Msg to PC
    ---------
    7 data bytes
        the first one is the box ID
        two are the firmware version (M, m)
        last 4 are the serial number

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_GET_CONFIGURATION
    STATE_IDLE
        Otherwise
    **/
    if (currentState != STATE_GET_CONFIGURATION){
        raise(ERROR_RUNTIME);
        return;
    }

    // Note that the ATmega32U4 of Arduino Micro uses
    // little-endian memory layout, i.e., LSB is
    // at lower memory index

	byte serial_nr[4];
	getSerialNR(serial_nr);
    byte configuration[7] = {BOX_ID,
                             FIRMWARE_VERSION_MAJOR,
                             FIRMWARE_VERSION_MINOR};
	int address = 0;
	while(address <= 3){
	  configuration[address + 3] = serial_nr[address];
	  address += 1;
	}
    encodeAndSend(configuration, LENGTH(configuration));
    currentState = STATE_IDLE;
}


/** Handler of STATE_SET_SERIAL_NR */
void setSerialNr(){
   /**
    Writes the assigned serial number to the EEPROM.

    Reads
    -----
    data_received

    Msg to PC
    ---------
    PC_OK

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_SERIAL_NR
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If the sent serial number is not 4 bytes long or if the
        sent data contains bytes that do not match the decimal
        ASCII representation of a capital letter or a number.
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_SERIAL_NR message
        and the receipt of data.
    STATE_SET_SERIAL_NR (stays)
        While waiting for data from the PC
    STATE_IDLE
        Successfully finished
    **/
    if (currentState != STATE_SET_SERIAL_NR){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage){  // waiting for data from the PC
        checkIfTimedOut();
        return;
    }

    // Data has arrived
    waitingForDataFromPC = false;
    newMessage = false;

    // Check that we got 4 bytes for the serial number.
    if (msgLength != 4){
        raise(ERROR_MSG_DATA_INVALID);
        return;
    }
	
	writeSerialNR(data_received);
	
    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}
