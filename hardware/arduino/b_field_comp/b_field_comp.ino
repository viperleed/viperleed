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

#include "states-def.h"                 // Basic state machine definitions
#include "arduino_utils.h"           // from '../lib'; for setChipSelectHigh
#include "viper-serial.h"             // Serial communication functions
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
    float current, temperature, shuntvoltage, voltage;                            // TODO: globals?
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

    compensateForSupplyVoltage();                                                    // TODO: Really here?
    compensateForThermalDrift();

    // Finally do what's needed in the currentState
    switch (currentState){
        case STATE_GET_CONFIGURATION:
            getConfiguration();
            break;
        case STATE_SET_SERIAL_NR:
            setSerialNr();
            break;
        case STATE_ERROR:
            handleErrors();
            break;
        case STATE_SET_DUTY_CYCLE:
            setDutyCycle();
            break;
        case STATE_SET_CURRENT:
            setCurrent();
            break;
        case STATE_MEASURE:
            measure();
            break;
        case STATE_SET_TIME_CONSTANT:
            setTimeConstant();
            break;
        case STATE_SET_ENABLE_DYNAMIC:
            setEnableDynamic();
            break;
        case STATE_SET_TRANSFORMATION_MATRIX:
            setTransformationMatrix();
            break;
        case STATE_SET_CALIBRATION_CURVE:
            setCalibrationCurve();
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
        case PC_RESET: break;
        case PC_SET_SERIAL_NR: break;
        case PC_STOP: break;
        case PC_SET_DUTY_CYCLE: break;
        case PC_SET_CURRENT: break;
        case PC_MEASURE: break;
        case PC_SET_TIME_CONSTANT: break;
        case PC_SET_ENABLE_DYNAMIC: break;
        case PC_SET_TRANSFORMATION_MATRIX: break;
        case PC_SET_CALIBRATION_CURVE: break;
       // case PC_: break;
        default:
            raise(ERROR_MSG_UNKNOWN);
            return false;
    }
    return true;                                                                // TODO: check that the received command is one of the commands that initiates one of the defined states
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
        case PC_RESET:
            encodeAndSend(PC_OK);
            reset();
            break;
        case PC_SET_SERIAL_NR:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_SERIAL_NR;
            break;
        case PC_STOP:
            currentState = STATE_IDLE;
            encodeAndSend(PC_OK);
            break;
        case PC_SET_DUTY_CYCLE:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_DUTY_CYCLE;
            break;
        case PC_SET_CURRENT:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_CURRENT;
            break;
        case PC_MEASURE:
            currentState = STATE_MEASURE;
            break;
        case PC_SET_TIME_CONSTANT:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_TIME_CONSTANT;
            break;
        case PC_SET_ENABLE_DYNAMIC:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_ENABLE_DYNAMIC;
            break;
        case PC_SET_TRANSFORMATION_MATRIX:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_TRANSFORMATION_MATRIX;
            break;
        case PC_SET_CALIBRATION_CURVE:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_CALIBRATION_CURVE;
            break;
        /*case PC_:
            currentState = STATE_;
            break;*/
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
    getSerialNumber(serial_nr);
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


/** Handler of STATE_ERROR */
void handleErrors(){
    /**Clean up after an error, and report it to the PC.                    //TODO: handle coils in case of error

    Reads
    -----
    currentState, errorTraceback

    Msg to PC
    ---------
    First message : PC_ERROR
    Second message : errorTraceback, i.e., two bytes
        First byte : state of the Arduino while the error occurred
        Second byte : error code

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within a STATE_ERROR
    STATE_IDLE
        Otherwise
    **/
    if (currentState != STATE_ERROR){
        raise(ERROR_RUNTIME);
        return;
    }

    // First, report the error, so the PC knows
    // there may be some cleanup going on
    encodeAndSend(PC_ERROR);
    encodeAndSend(errorTraceback, LENGTH(errorTraceback));

    // Then clean up possible mess that caused the error
    switch(errorTraceback[1]){
        case ERROR_SERIAL_OVERFLOW:
            // Discard all characters in the serial buffer,
            // since messages are anyway likely corrupt.
            while(Serial.available()) Serial.read();
            break;
        case ERROR_MSG_TOO_LONG:
            numBytesRead = 0;
            readingFromSerial = false;
            break;
    }
    waitingForDataFromPC = false;
    newMessage = false;
    currentState = STATE_IDLE;
}


/** Handler of PC_RESET */
void reset(){                                      //TODO: add reset settings for coils
    /**
    Reset the Arduino, mimicking as much as possible the state right
    after boot-up or a hardware reset.

    Writes
    ------
    numBytesRead, msgLength, readingFromSerial, newMessage,
    waitingForDataFromPC,

    Goes to state
    -------------
    STATE_IDLE
    **/
    numBytesRead = 0;
    msgLength = 0;
    readingFromSerial = false;
    newMessage = false;

    currentState = STATE_IDLE;
    waitingForDataFromPC = false;

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

    storeSerialNumber(data_received);

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}


/** Handler of STATE_SET_DUTY_CYCLE */
void setDutyCycle(){
     /**
    Wait until we get data from the PC and set a duty cycle for the respective coil.

    For each coil, we need 4 bytes from the data_received[] buffer.
    Their meaning is:
    - Set duty cycle for coil 1
    - Set duty cycle for coil 2

    Reads
    -----
    currentState, newMessage, waitingForDataFromPC, data_received, msgLength

    Writes
    ------
    newMessage, waitingForDataFromPC, dutyCycle,

    Msg to PC
    ---------
    PC_OK, signaling that we're done setting the duty cycle

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_DUTY_CYCLE
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_DUTY_CYCLE message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_SET_DUTY_CYCLE){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // We need exactly 8 bytes, 4 bytes for each float representing one duty cycle
        if (msgLength != 8){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }
        for(int j=0; j<2;j++){
            for(int i = 0; i < 4 ;i++){
                dutyCycle[j].asBytes[4-i] = data_received[i+j*4];
             }
        }
                                                                      //TODO: catch errors from coils and report them
        coil_1.set_duty_cycle(dutyCycle[0].asFloat);     // This sets the duty cycle for coil 1
        coil_2.set_duty_cycle(dutyCycle[1].asFloat);     // This sets the duty cycle for coil 2

        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
     }
}


/** Handler of STATE_SET_CURRENT */
void setCurrent(){
   /** Wait until we get data (a target current) from the PC, convert the received data in order to set a duty cycle 

    For each coil, we need 4 bytes from the data_received[] buffer.
    Their meaning is:
    - Target current for coil 1
    - Target current for coil 2

    Reads
    -----
    currentState, newMessage, waitingForDataFromPC, data_received, msgLength

    Writes
    ------
    newMessage, waitingForDataFromPC, targetCurrent,

    Msg to PC
    ---------
    PC_OK, signaling that we're done setting the current/duty cycle

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_CURRENT
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_CURRENT message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_SET_CURRENT){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // 2 floats -> msgLength is 8
        if (msgLength != 8){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }
        for(int j=0; j<2;j++){
            for(int i = 0; i < 4 ;i++){
                targetCurrent[j].asBytes[3-i] = data_received[i+j*4];
             }
        }

        //Function convertCurrentToDutyCycle() needs to be written
        convertCurrentToDutyCycle();

        coil_1.set_duty_cycle(dutyCycle[0].asFloat);     // This sets the duty cycle for coil 1 using the current
        coil_2.set_duty_cycle(dutyCycle[1].asFloat);     // This sets the duty cycle for coil 2 using the current

        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
     }
}


void convertCurrentToDutyCycle(){
    //This function converts the target current given by the PC, to a duty cycle which will be set on the respective coils
    float x = 0;
    dutyCycle[0].asFloat = targetCurrent[0].asFloat * x;        // TODO: implement
    dutyCycle[1].asFloat = targetCurrent[1].asFloat * x;        // TODO: implement
}


/** Handler of STATE_MEASURE */
void measure(){
   /** We are measuring the current and the supply voltage.

    Reads
    -----
    currentState

    Writes
    ------
    voltage1, voltage2, current1, current2

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_MEASURE
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_MEASURE){
        raise(ERROR_RUNTIME);
        return;
    }

    floatOrBytes voltage1;  //TODO: make implementation of the returning data more efficient if this works
    floatOrBytes voltage2;
    floatOrBytes current1;
    floatOrBytes current2;

    //Measure the current of the respective coil
    current1.asFloat = coil_1.get_current();
    current2.asFloat = coil_2.get_current();

    //Measure the voltage at pin8 (motor driver's supply voltage) of the respective coil
    voltage1.asFloat = coil_1.get_voltage();
    voltage2.asFloat = coil_2.get_voltage();


    sendFloatToPC(current1);
    sendFloatToPC(current2);
    sendFloatToPC(voltage1);
    sendFloatToPC(voltage2);



    currentState = STATE_IDLE;
}


void sendFloatToPC(floatOrBytes value){
    // Since the ATMega32u4 (Arduino Micro) uses little-endian memory layout, we
    // have flip the bytes over to maintain consistency of our messages, which
    // are in big-endian order
    byte littleToBigEndian[4];

    for (int i = 0; i < 4; i++){
        littleToBigEndian[i] = value.asBytes[3-i];
    }
    encodeAndSend(littleToBigEndian, 4);
}


/** Handler of STATE_SET_TIME_CONSTANT */
void setTimeConstant(){
   /** Wait until we get data (time constant) from the PC, convert the received data

    For each coil, we need 4 bytes from the data_received[] buffer. ???
    Their meaning is:
    - ...
    - ...

    Reads
    -----
    currentState, newMessage, waitingForDataFromPC, data_received, msgLength

    Writes
    ------
    newMessage, waitingForDataFromPC, timeConstant

    Msg to PC
    ---------
    PC_OK, signaling that we're done

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_TIME_CONSTANT
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_TIME_CONSTANT message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_SET_TIME_CONSTANT){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // 2 floats ->  msgLength is 8
        if (msgLength != 8){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }

        for(int j=0; j<2;j++){
            for(int i = 0; i < 4 ;i++){
                timeConstant[j].asBytes[3-i] = data_received[i+j*4];
             }
        }

        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
     }
}


/** Handler of STATE_SET_ENABLE_DYNAMIC */
void setEnableDynamic(){
   /** Wait until we get data from the PC, convert the received data

    For each coil, we need 4 bytes from the data_received[] buffer.  ???
    Their meaning is:
    - ...
    - ...

    Reads
    -----
    currentState, newMessage, waitingForDataFromPC, data_received, msgLength

    Writes
    ------
    newMessage, waitingForDataFromPC, enableDynamic

    Msg to PC
    ---------
    PC_OK, signaling that we're done

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_ENABLE_DYNAMIC
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_ENABLE_DYNAMIC message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_SET_ENABLE_DYNAMIC){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // Our msgLength is 1 as we get one byte.
        if (msgLength != 1){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }

        if (data_received[0] == 1){
            enableDynamic = true;
        }

         if (data_received[0] == 0){
            enableDynamic = false;
        }

        else{
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }

        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
     }
}


/** Handler of STATE_SET_TRANSFORMATION_MATRIX */
void setTransformationMatrix(){
   /** Wait until we get data from the PC (what data?), convert the received data in order to set a transformation matrix

    For each coil, we need ? bytes from the data_received[] buffer.
    Their meaning is:
    -
    -

    Reads
    -----
    currentState, newMessage, waitingForDataFromPC, data_received, msgLength

    Writes
    ------
    newMessage, waitingForDataFromPC, transformationMatrix

    Msg to PC
    ---------
    PC_OK, signaling that we're done

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_TRANSFORMATION_MATRIX
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_TRANSFORMATION_MATRIX message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_SET_TRANSFORMATION_MATRIX){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // 6 floats in the matrix -> msgLength is 24
        if (msgLength != 24){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }
        for(int m=0; m<3;m++){
            for(int j=0; j<2;j++){
                for(int i = 0; i < 4 ;i++){
                    transformationMatrix[m][j].asBytes[3-i] = data_received[i+j*4+m*8];
                 }
            }
        }
        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
     }
}


/** Handler of STATE_SET_CALIBRATION_CURVE */
void setCalibrationCurve(){
   /** Wait until we get data from the PC (what data?), convert the received data in order to set a calibration curve

    For each coil, we need ? bytes from the data_received[] buffer.
    Their meaning is:
    -
    -

    Reads
    -----
    currentState, newMessage, waitingForDataFromPC, data_received, msgLength

    Writes
    ------
    newMessage, waitingForDataFromPC, calibrationCurve

    Msg to PC
    ---------
    PC_OK, signaling that we're done

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_CALIBRATION_CURVE
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_CALIBRATION_CURVE message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_IDLE  :
        Successfully finished

    **/

     if (currentState != STATE_SET_CALIBRATION_CURVE){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // N bytes in the array -> msgLength is N
        if (msgLength != N){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }

        for(int j=0; j<N;j++){
            calibrationCurve[j] = data_received[j];
        }

        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
     }
}