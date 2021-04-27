/*
ViPErLEED - Firmware for Arduino hardware controller
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 16.04.2021
---------------------
*/

// Libraries
#include <Arduino.h>
#include <SPI.h>

#include "viper-ino.h"   // Arduino-related settings. Includes ADC and DAC

#define DEBUG              true      // Debug mode, writes to serial line, for use in serial monitor

// Firmware version (MAX: v255.255). CURENTLY: v0.1
#define FIRMWARE_VERSION_MAJOR    0  // max 255
#define FIRMWARE_VERSION_MINOR    1  // max 255






/** ---------------------------- INITIALIZATION ---------------------------- **/

void setup() {
    /** Set up the Arduino board ans all the hardware.

    The setup routine runs once on power-up or on hardware reset.
    */
    #ifdef DEBUG
        // initialize serial communication (emulation
        // on USB) at 9600 bits per second for debug
        delay(1000);
        Serial.begin(9600);
        delay(1000);
    #endif

    analogReference(INTERNAL);  // Used for LM35 and for checking relay presence

    // Set all inverted chip-select pins to high via pull-up (to avoid low glitch)
    setChipSelectHigh(CS_DAC);
    setChipSelectHigh(CS_ADC_0);
    setChipSelectHigh(CS_ADC_1);
    pinMode(SCK, OUTPUT);  // Should not be needed, but it did not work without
    pinMode(MOSI, OUTPUT);

    // Reset ADC and DAC, set DAC output to zero,
    // and set variables to well-defined defaults
    reset();
}






/** ----------------------------- STATE MACHINE ---------------------------- **/

void loop() {
    /**Read from serial line and handle the state machine.**/
    // See if there is anything coming from the PC.
    readFromSerial();

    // Then see if the current state needs updating,
    // depending on the message received.
    updateState();

    // Finally do what's needed in the currentState
    switch (currentState){
        case STATE_GET_CONFIGURATION:
            getConfiguration();
            break;
        case STATE_CALIBRATE_ADCS:
            calibrateADCsAtAllGains();
            break;
        case STATE_SET_UP_ADCS:
            prepareADCsForMeasurement();
            break;
        case STATE_SET_VOLTAGE:
            setVoltage();
            break;
        case STATE_TRIGGER_ADCS:
            waitAndTriggerMeasurements();
            break;
        case STATE_MEASURE_ADCS:
            measureADCs();
            break;
        case STATE_ADC_VALUES_READY:
            sendMeasuredValues();
            break;
        case STATE_AUTOGAIN_ADCS:
            findOptimalADCGains();
            break;
        case STATE_ERROR:
            handleErrors();
            break;
    }
}


void updateState() {
    /**
    If there is an unprocessed serial message, decide whether
    this requires us to change the state of the Arduino, and
    do some preparation that may be needed before entering
    the new state.
    
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
    if (not newMessage) return;

    switch(data_received[0]){
        case PC_CONFIGURATION:
            initialTime = millis();
            currentState = STATE_GET_CONFIGURATION;
            break;
        case PC_CALIBRATION:
            waitingForDataFromPC = true;
            initialTime = millis();
            calibrationGain = 0;
            currentState = STATE_CALIBRATE_ADCS;
            break;
        case PC_SET_UP_ADCS:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_UP_ADCS;
            break;
        case PC_SET_VOLTAGE:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_VOLTAGE;
            break;
        case PC_AUTOGAIN:
            initialTime = millis();
            prepareForAutogain();
            currentState = STATE_AUTOGAIN_ADCS;
            break;
        case PC_TRIGGER_ADCS:
            initialTime = millis();
            currentState = STATE_TRIGGER_ADCS;
            break;
        case PC_RESET:
            reset();
            break;
    }
    newMessage = false;
}


bool didTimeOut(){
    /**Return whether the Arduino waiting has timed out.
    
    Reads
    -----
    initialTime
    
    Returns
    -------
    true if timed out
    
    Goes to state
    -------------
    STATE_ERROR with ERROR_TIMEOUT
        If timed out
    **/
    if((millis() -  initialTime) > TIMEOUT){
        errorTraceback[0] = currentState;
        errorTraceback[1] = ERROR_TIMEOUT;
        currentState = STATE_ERROR;
        return true;
        }
    return false;
}


bool notInStateAndError(byte state){
    /**Return true if we are not in the given state.
    
    Also brings the system to a RUNTIME_ERROR.
    
    Reads
    -----
    currentState
    
    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If not in state
    **/
    if (currentState != state){
        errorTraceback[0] = currentState;
        errorTraceback[1] = ERROR_RUNTIME;
        currentState = STATE_ERROR;
        return true;
    }
    return false;
}






/** ----------------------------- COMMUNICATION ---------------------------- **/

void readFromSerial() {
    /**
    Store bytes (i.e., characters) read from the serial line (i.e., PC)
    in serialInputBuffer. A message is considered complete when we
    receive a MSG_END. When this happens, serialInputBuffer will be
    [MSG_START, byte with length of decoded message, message, MSG_END]
    
    A full message will likely be read during a single call to this
    function, unless its characters come in very slowly. In this case
    it will take a few loop iterations to read a full message.

    Writes
    ------
    numBytesRead, serialInputBuffer, readingFromSerial, msgLength

    Goes to state
    -------------
    STATE_ERROR + ERROR_SERIAL_OVERFLOW
        In case the Arduino serial buffer reaches its limit
    STATE_ERROR + ERROR_MSG_INCONSITENT
        In case the number of bytes read does not fit with the
        number expected from the first byte after MSG_START

    */
    // Do something only if there is data on the serial line
    if(not Serial.available())  return;

    if(Serial.available() >= SERIAL_BUFFER_SIZE){  // Should never be '>', but better safe than sorry
        // The serial buffer is full and it potentially
        // already discarded some of the bytes that came
        errorTraceback[0] = currentState;
        errorTraceback[1] = ERROR_SERIAL_OVERFLOW;
        currentState = STATE_ERROR;
        // The buffer will be flushed in the error handler 
        return;
    }

    while (Serial.available()){
        byte byteRead = Serial.read();

        // New message
        if (byteRead == MSG_START) {
            numBytesRead = 0;
            readingFromSerial = true;
        }

        // Accumulate characters
        if(readingFromSerial) {
            // Make sure we are not going to write
            // past the end of serialInputBuffer
            if (numBytesRead == MSG_MAX_LENGTH) {
                errorTraceback[0] = currentState;
                errorTraceback[1] = ERROR_MSG_TOO_LONG;
                currentState = STATE_ERROR;
                return;
            }
            serialInputBuffer[numBytesRead] = byteRead;
            numBytesRead++;
        }

        // A full message has been read
        if (byteRead == MSG_END) {
            readingFromSerial = false;
            newMessage = decodeAndCheckMessage();
            break;
        }
        
        // Delay a very little bit to make sure new characters
        // come in, if any. This should be enough to read in a
        // whole message. If it isn't the case, we will finish
        // during the next state-loop iteration anyway.
        delayMicroseconds(20);
    }
}


bool decodeAndCheckMessage(){
    /**
    Move the interesting bytes of serialInputBuffer[]
    into data_received[], and return whether the
    message read is acceptable.

    In practice, only the actual characters are kept.
    This means:
    (1) Skipping:
           MSG_START == serialInputBuffer[0]
           no. of bytes after decoding == serialInputBuffer[1]
           MSG_END == serialInputBuffer[last]
    (2) Decoding bytes with value MSG_SPECIAL_BYTE. In this case, the
        actual character is MSG_SPECIAL_BYTE + the next character.

    A message is considered acceptable if:
    (1) the length contained in the message fits the number
        of bytes decoded
    (2) when it is a one-character message (i.e., a command)
        it should be one of the known commands
    (3) when it is longer (i.e., it's data), we should be
        in a state that expects data

    If a message is not acceptable, it also brings the system
    to a STATE_ERROR with the appropriate error code

    Reads
    -----
    serialInputBuffer

    Writes
    ------
    data_received, currentState, msgLength

    Returns
    -------
    True if the message is acceptable

    Goes to state
    -------------
    (unchanged)
        if message is acceptable
    STATE_ERROR with ERROR_MSG_INCONSITENT
        if the number of decoded bytes does not match the
        length expected from the value that came with the
        message itself
    STATE_ERROR with ERROR_MSG_UNKNOWN
        if the message is an unknown command, or if we
        got some 'data' while we were not expecting any
    */
    // Decode the message, starting at the second
    // byte, and going up to numBytesRead-1. This
    // skips MSG_START, the length, and MSG_END
    byte numDecodedBytes = 0;
    for (byte nthByte = 2; nthByte < numBytesRead; nthByte++) {
        byte decodedByte = serialInputBuffer[nthByte];
        if (decodedByte == MSG_SPECIAL_BYTE) {
            // The actual character is MSG_SPECIAL_BYTE + the next byte
            nthByte++;
            decodedByte += serialInputBuffer[nthByte];
        }
        data_received[numDecodedBytes] = decodedByte;
        numDecodedBytes++;
    }

    // Check that the number of bytes decoded fits
    msgLength = serialInputBuffer[1];
    if (msgLength != numDecodedBytes){
        errorTraceback[0] = currentState;
        errorTraceback[1] = ERROR_MSG_INCONSITENT;
        currentState = STATE_ERROR;
        return false;
        }

    if (numDecodedBytes > 1){
        // Message is some data
        if (not waitingForDataFromPC){
            // But we're not expecting any
            errorTraceback[0] = currentState;
            errorTraceback[1] = ERROR_MSG_UNKNOWN;
            currentState = STATE_ERROR;
            return false;
        }
        // Defer checking to state handlers
        return true;
    }

    // Check that it is one of the understandable commands
    switch(data_received[0]){
        case PC_AUTOGAIN: break;
        case PC_CALIBRATION: break;
        case PC_CONFIGURATION: break;
        case PC_SET_UP_ADCS: break;
        case PC_TRIGGER_ADCS: break;
        case PC_RESET: break;
        case PC_SET_VOLTAGE: break;
        default:
            errorTraceback[0] = currentState;
            errorTraceback[1] = ERROR_MSG_UNKNOWN;
            currentState = STATE_ERROR;
            return false;
    }
    return true;
}


void encodeAndSend(byte singleByte){
/*
 * Prepares message before sending it to the PC. Puts single
 * byte into an array and forwards the array to the "real"
 * encode message
 * This overloaded function essentially prepares a single-byte-long message
 * to be actually encoded in the next function. Having the two with the same
 * names prevents the rest of the code from having to figure out which function
 * to call depending on whether the message is a single byte or a byte array
 *
 * Parameters:
 * -----------
 * singleByte : single byte
 *
 * Return into function:
 * -----------
 * NewVar[]: byte array with one value
 */
  byte byteArray[1];
  byteArray[0] = singleByte;
  encodeAndSend(byteArray, 1);
}


void encodeAndSend(byte *byteArray, int len){
/*
 * Prepares message before sending it to the PC. Changes every
 * byte which happens to have the same value as a MSG_START, an MSG_END or
 * a MSG_SPECIAL_BYTE to two bytes with a leading "MSG_SPECIAL_BYTE" and a following
 * "byte - MSG_SPECIAL_BYTE."
 * The message is put into the array "data_send[]". The length
 * of the array is defined in the variable "numBytesBeforeEncoding".
 *
 * Parameters:
 * -----------
 * byteArray : byte array
 * len : integer, length of byte array
 *
 * Returns into globals
 * -----------
 * data_send: byte array
 */
  byte numBytesAfterEncoding = 0;
  byte numBytesBeforeEncoding = 0; // !!!!!!! Can we remove this? Byte cast not clear
  for(int i = 0; i < len; i ++){
    if(byteArray[i] >= MSG_SPECIAL_BYTE){
      data_send[numBytesAfterEncoding] = MSG_SPECIAL_BYTE;
      numBytesAfterEncoding++;
      data_send[numBytesAfterEncoding] = byteArray[i] - MSG_SPECIAL_BYTE;
    }
    else{
      data_send[numBytesAfterEncoding] = byteArray[i];
    }
    numBytesAfterEncoding++;
  }
  numBytesBeforeEncoding = len;
/*
 * Sends byte array "data_send" (i.e., the actual message) to PC as:
 *   MSG_START
 *   numbers of bytes in actual message
 *   actual message
 *   MSG_END
 */
  Serial.write(MSG_START);
  Serial.write(numBytesBeforeEncoding);
  Serial.write(data_send, numBytesAfterEncoding);
  Serial.write(MSG_END);
}






/** -------------------- PREPARE FOR ENTERING NEW STATES ------------------- **/

void prepareForAutogain(){
    /**Prepare the system to enter a STATE_AUTOGAIN_ADCS.
    
    Writes
    ------
    adc0Gain, adc1Gain, numMeasurementsDone, summedMeasurements,
    numMeasurementsToDoBackup, numMeasurementsToDo
    **/
    // Reset data, as we will store measurements in the same arrays
    resetMeasurementData();
    
    // Gains to zero, where we will measure for the gain optimization
    adc0Gain = 0;
    adc1Gain = 0;
    
    // Calibrate at the high frequency used for autogain
    selfCalibrateAllADCs(AD7705_500HZ);  // Takes ~13 ms
    
    // Backup the number of measurement that will be
    // used during normal measurements later on. This
    // is restored at the end of the gain optimization
    numMeasurementsToDoBackup = numMeasurementsToDo;
    numMeasurementsToDo = 25;
}






/** ---------------------- STATE & PC-REQUEST HANDLERS --------------------- **/

/** Handler of STATE_GET_CONFIGURATION */
void getConfiguration(){
    /**Send firmware version and hardware configuration to PC.

    Writes
    -----
    hardwareDetected

    Msg to PC
    ---------
    6 data bytes
        first four are firmware version (M, M, m, m)
        last two are hardware configuration as bitmask

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_GET_CONFIGURATION
    STATE_IDLE
        Otherwise
    **/
    if (notInStateAndError(STATE_GET_CONFIGURATION))
        return;

    hardwareDetected.asInt = getHardwarePresent();
    byte configuration[4] = {FIRMWARE_VERSION_MAJOR,
                             FIRMWARE_VERSION_MINOR,
                             hardwareDetected.asBytes[0],
                             hardwareDetected.asBytes[1]};
    encodeAndSend(configuration, 6);
    currentState = STATE_IDLE;
}


// For reasons unknown, the Arduino compiler complains if this
// function is declared after its caller (the next one). Having
// a prototype in the .h file does not help.
void storeAllSelfCalibrationResults(int32_t targetArray[N_MAX_ADCS_ON_PCB][2]) {
    /**
    Store the result of the last self calibration
    of both ADCs (if present) to the array passed

    Parameters
    ----------
    targetArray : 2x2 int32_t
        The array where the calibration results are to be saved.
        The first index identifies which ADC, the second one is
        used to store offset (at index 0) and gain (at index 1).

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel
    **/
    for (int iADC = 0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
        // bool adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;       // BUG: ADC_0_PRESENT/ADC_1_PRESENT are uint16_t, not bool!
        uint16_t adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;      // TODO: using the ADCs struct would make it easier, as you would just need to do externalADCs[iADC].present
        if (hardwareDetected.asInt & adcPresentMask) {
            byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;                    // TODO: similarly, with the ADCs struct, would be enough to do externalADCs[iADC].chipSelect
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: same thing here: externalADCs[iADC].channel
            int32_t offs = AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET);
            int32_t gain = AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN);
            targetArray[iADC][0] = offs;
            targetArray[iADC][1] = gain;
        }
    }
}


/** Handler of STATE_CALIBRATE_ADCS */
void calibrateADCsAtAllGains(){
    /**Calibrate (in parallel) the available ADCs for all possible gains.

    The calibration results are stored, so they can be fetched
    later to make gain-switching faster.
    
    Calibration is performed for those channels and update
    rate that come from the PC as the first data message
    after PC_CALIBRATION.

    Each gain value is done over a separate state-loop. This means
    that it takes approximately 360 ms to finish a single run of
    this function (120 ms x 3 points for computing medians), i.e.,
    that the serial line remains unresponsive for ~360 ms.

    Should one need to change the updateRate used during measurements,
    the system must run again through this STATE_CALIBRATE_ADCS.

    Reads
    -----
    msgLength

    Writes
    ------
    initialTime, newMessage, waitingForDataFromPC,
    adc0Channel, adc1Channel, adcUpdateRate,
    calibrationGain, adc0Gain, adc1Gain, selfCalDataVsGain

    Msg to PC
    ---------
    PC_OK.

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_CALIBRATE_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If the ADC channels and update rate data needed from
        the PC do not arrive within 5 s of entering the state
    STATE_ERROR : ERROR_TIMEOUT
        If calibration for a single state-loop iteration (i.e.,
        one value of gain) takes longer than 5 s
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we get the wrong number of data bytes, or if either
        one of the channels or the update rate is not one of
        the acceptable values
    STATE_CALIBRATE_ADCS (stays)
        Until all the gain values have been calibrated
        for the currently selected channels
    STATE_IDLE
        After calibration is done
    */
    if (notInStateAndError(STATE_CALIBRATE_ADCS))
        return;

    if(not newMessage){  // waiting for data from the PC
        if (didTimeOut())
            waitingForDataFromPC = false;
        return;
    }
    else {
        // Data has arrived. This part runs only once,
        // the first time the data arrives

        // Check that we got the data that we expected, i.e.,
        // (1) 2 bytes for channels, one for update rate
        if (msgLength != 3){
            errorTraceback[0] = currentState;
            errorTraceback[1] = ERROR_MSG_DATA_INVALID;
            currentState = STATE_ERROR;
            waitingForDataFromPC = false;
            newMessage = false;
            return;
        }
        
        // (2) The channels are acceptable
        for (byte i = 0; i < 2; i++){
            byte channel = data_received[i];
            if (channel != AD7705_CH0 or channel != AD7705_CH1){
                errorTraceback[0] = currentState;
                errorTraceback[1] = ERROR_MSG_DATA_INVALID;
                currentState = STATE_ERROR;
                waitingForDataFromPC = false;
                newMessage = false;
                return;
            }
        }
        
        adc0Channel = data_received[0];
        adc1Channel = data_received[1];
        
        // (2) The update rate is acceptable
        if (not isAcceptableADCUpdateRate(data_received[2])){
            errorTraceback[0] = currentState;
            errorTraceback[1] = ERROR_MSG_DATA_INVALID;
            currentState = STATE_ERROR;
            waitingForDataFromPC = false;
            newMessage = false;
            return;
        }
        adcUpdateRate = data_received[2];
        
    }

    // Data is OK
    newMessage = false;
    waitingForDataFromPC = false;
    
    // May timeout in the following if the ADC power is disconnected
    // while we're waiting for a calibration to be finished
    initialTime = millis();

    // The next part runs once every state loop
    if (calibrationGain <= AD7705_MAX_GAIN) {        
        // Prepare a place to store values from which medians are
        // calculated: three for median, N_MAX_ADCS_ON_PCB for ADCs,
        // last index is offset(0) & gain(1)
        int32_t selfCalDataForMedian[3][N_MAX_ADCS_ON_PCB][2];

        if (didTimeOut()) {
            adc0Gain = 0;
            adc1Gain = 0;
            return;
        }

        adc0Gain = calibrationGain;
        adc1Gain = calibrationGain;
        // Calibrate 3 times each of the ADCs (in parallel) at
        // the current calibrationGain, and for the adc0Channel
        // and adc1Channel currently selected.
        for (int i = 0; i < 3; i++) {
            selfCalibrateAllADCs(adcUpdateRate);
            storeAllSelfCalibrationResults(selfCalDataForMedian[i]);
            #ifdef DEBUG
                Serial.print("gain=");
                Serial.print(calibrationGain);
                Serial.print(" cOffs=");
                Serial.print(selfCalDataForMedian[i][0][0]);
                Serial.print(" cGain=");
                Serial.println(selfCalDataForMedian[i][0][1]);
            #endif
        }

        // Compute the median of the values above, and store
        // it for each ADC, so it can be fetched later
        for (int iADC=0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: using the externalADCs array of struct would make this easier: externalADCs[iADC].channel
            for (int offsetOrGain = 0; offsetOrGain < 2; offsetOrGain++) {
                int32_t median = getMedian32(
                    selfCalDataForMedian[0][iADC][offsetOrGain],
                    selfCalDataForMedian[1][iADC][offsetOrGain],
                    selfCalDataForMedian[2][iADC][offsetOrGain]
                    );
                selfCalDataVsGain[calibrationGain][iADC][channel][offsetOrGain] = median;
            }
        }
        calibrationGain++;
    }

    if(calibrationGain > AD7705_MAX_GAIN){
        // Calibration is over.
        // Remember which channels we have calibrated
        for (int iADC=0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
            byte channel = iADC==0 ? adc0Channel : adc1Channel;
            calibratedChannels[iADC][channel] = true;
        }

        // Set gains to the lowest value, and load the
        // new gains and their calibration in the ADCs
        adc0Gain = 0;
        adc1Gain = 0;
        setAllADCgainsAndCalibration();

        // And tell the PC that we're done
        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
    }
}


/** Handler of STATE_SET_UP_ADCS */
void prepareADCsForMeasurement(){
    /**
    Initialize the ADC from the parameters stored in the first 5 bytes of
    data_received. The bytes have the following meaning:
    - [0] channel of ADC#0
    - [1] channel of ADC#1
    - [2 + 3] number of measurements that will be averaged. [2] is the
              high byte, [3] the low byte. Data is interpreted as a uint16_t.

    - // TODO: remove from Python: measurement frequency for ADCs #0 and #1
    - // TODO: remove from Python: maximum gain value for ADC #0 and #1

    Reads
    -----
    data_received

    Writes
    ------
    newMessage, adc0Channel, adc1Channel, adcUpdateRate, numMeasurementsToDo,
    currentState

    Msg to PC
    ---------
    PC_OK

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_UP_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_UP_ADCS message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If either channel is an acceptable values
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_SET_UP_ADCS (stays)
        While waiting for data from the PC
    STATE_IDLE
        Successfully finished
    **/
    if (notInStateAndError(STATE_SET_UP_ADCS))
        return;

    if (not newMessage){  // waiting for data from the PC
        if(didTimeOut()){
            waitingForDataFromPC = false;
            newMessage = false;
        }
        return;
    }
    
    // Data has arrived
    waitingForDataFromPC = false;
    newMessage = false;

    // Do some checking on the data we received:
    // (1) Correct number of bytes? [2(channels) + 2(no. measurements)]
    if (msgLength != 4){
        errorTraceback[0] = currentState;
        errorTraceback[1] = ERROR_MSG_DATA_INVALID;
        currentState = STATE_ERROR;
        return;
    }
    
    // (2) channels acceptable?
    for (byte i = 0; i < 2; i++){
        byte channel = data_received[i];
        if (channel != AD7705_CH0 or channel != AD7705_CH1){
            errorTraceback[0] = currentState;
            errorTraceback[1] = ERROR_MSG_DATA_INVALID;
            currentState = STATE_ERROR;
            return;
        }
    }

    adc0Channel = data_received[0];
    adc1Channel = data_received[1];

    numMeasurementsToDo = data_received[3] << 8 | data_received[4];

    setAllADCgainsAndCalibration();
    if (currentState == STATE_ERROR)   // Some channel was not calibrated
            return;

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}


/** Handler of STATE_SET_VOLTAGE */
void setVoltage(){
    /**
    Ask the DAC to provide a new voltage, if new settings are available.

    Needs 4 bytes from the data_received[] buffer. The meaning is:
    - [0(MSB) + 1(LSB)] DAC value to be set (uint16_t).
    - [2(MSB) + 3(LSB)] dacSettlingTime (uint16_t). Time interval to
            wait before considering the DAC value stable

    Reads
    -----
    data_received

    Writes
    ------
    dacSettlingTime, currentState, newMessage, summedMeasurements,
    numMeasurementsDone, adc0Gain, adc1Gain, adc0ShouldDecreaseGain,
    adc1ShouldDecreaseGain

    Msg to PC
    ---------
    None.

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_VOLTAGE
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_VOLTAGE message
        and the receipt of data
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_SET_VOLTAGE (stays)
        While waiting for data from the PC
    STATE_TRIGGER_ADCS
        Successfully finished
    **/
    if (notInStateAndError(STATE_SET_VOLTAGE))
        return;

    if(not newMessage){   // Waiting for data from PC
        if(didTimeOut())
            waitingForDataFromPC = false;
        return;
    }
    
    waitingForDataFromPC = false;

    uint16_t dacValue = data_received[0] << 8 | data_received[1];
    dacSettlingTime = data_received[2] << 8 | data_received[3];
    AD5683setVoltage(CS_DAC, dacValue);
    currentState = STATE_TRIGGER_ADCS;
    initialTime = millis();
    newMessage = false;
    resetMeasurementData();
    decreaseADCGainsIfNeeded();
}


/** Handler of STATE_TRIGGER_ADCS */
void waitAndTriggerMeasurements() {
    /**
    Wait for the voltage output to be stable, then
    trigger the (available) ADCs to start converting

    Reads
    -----
    hardwareDetected, adc0Channel, adc0Gain, adc1Channel, adc1Gain

    Writes
    ------
    currentState

    Msg to PC
    ---------
    PC_OK, signaling that the we will then begin collecting measurements

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_TRIGGER_ADCS
    STATE_TRIGGER_ADCS (stays)
        until the voltage output can be considered stable
    STATE_MEASURE_ADCS
        after the voltage output can be considered stable
    **/
    if (notInStateAndError(STATE_TRIGGER_ADCS))
        return;
    
    // Wait in the same state till the voltage output is stable
    if((millis() - initialTime) < dacSettlingTime) return;
    
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);
    
    // Signal the PC that we are now going to start the measurements
    encodeAndSend(PC_OK);
    
    // Switch over to measuring state
    currentState = STATE_MEASURE_ADCS;
}


/** Handler of STATE_MEASURE_ADCS */
void measureADCs(){
    /**Acquire measurements from the ADCs.

    (1) Take a defined number of measurements (as per numMeasurementsToDo)
        that are added together in summedMeasurements[].
    (2) When the number of measurement points reaches the requested amount,
        the state of the Arduino is changed to STATE_ADC_VALUES_READY, which
        will handle the communication of the measurement to the PC.

    Reads
    -----
    numMeasurementsDone, numMeasurementsToDo, initialTime

    Writes
    ------
    summedMeasurements, currentState, newMessage

    Msg to PC
    ---------
    None.

    Goes to state
    -------------
    STATE_ERROR : ERROR_TIMEOUT
        If it takes longer than 5s between the beginning of the state
        and completing the number of measurements that need to be averaged
    STATE_ERROR : ERROR_ADC_SATURATED
        If one of the inputs of the ADCs has reached saturation, and
        there is no room to decrease the gain.
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_MEASURE_ADCS (stays)
        Until finished measuring all the data points
    STATE_ADC_VALUES_READY
        Successfully finished
    **/
    if (notInStateAndError(STATE_MEASURE_ADCS))
        return;

    makeAndSumMeasurements();
    // TODO: if, by chance, the last value measured by one of the ADCs
    //       caused a solid saturation, this is not communicated in any
    //       way to the PC, since the rest of the code goes through, and
    //       the Arduino goes to STATE_ADC_VALUES_READY, swallowing
    //       the error. It could be solved easily by a simple check:
    //       if (currentState == STATE_ERROR) return;
    if(numMeasurementsDone == numMeasurementsToDo){
        currentState = STATE_ADC_VALUES_READY;
        return;
    }

    if(didTimeOut())
        newMessage = false;
}


/** Handler of STATE_ADC_VALUES_READY */
void sendMeasuredValues(){
    /**Send measurements back to the PC.

    Currently, it actually sends THREE separate messages. Each of
    the messages contains 4 data bytes that, after being re-packed
    (MSB first) correspond to the analog values measured in physical
    units like so:
    - I0(@LEED)   ADC#0, ch0: Volts
    - HV          ADC#0, ch1: Volts
    - I0(@sample) ADC#1, ch0: uAmps
    - AUX         ADC#1, ch1: Volts
    - LM35 (if present)     : degC
    - LM35 (if not present) : Volts

    The sequence of returned results is: LM35, ADC#0, ADC#1

    Notice that the I0 value (measured at the LEED electronics, i.e.,
    ADC#0, channel 0) will be returned in Volts, thus, the conversion
    to uAmps is delegated to the PC.

    Reads
    -----
    fDataOutput
    (indirectly) hardwareDetected, summedMeasurements, voltsPerBit,
                 numMeasurementsToDo, adc0Channel, adc1Channel,
                 adc0Gain, adc1Gain

    Writes
    ------
    (indirectly) fDataOutput, summedMeasurements, numMeasurementsDone

    Msg to PC
    ---------
    Three serial messages, each contains 4 bytes (+encoding) with
    values from LM35, ADC#0, ADC#1 (in this order)

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_ADC_VALUES_READY
    STATE_IDLE
        Otherwise
    **/
    // TODO: We may later remove this, if we want the option
    // of sending incomplete data back upon request from the PC
    if (notInStateAndError(STATE_ADC_VALUES_READY))
        return;
    
    // TODO: Ideally, one would like to rather return a SINGLE MESSAGE
    //       containing 3*4 significant bytes (+ encoding).
    //       In this case, the worst-case scenario message length would be
    //       N_MAX_MEAS*4(bytes)*2(encoding) + 3 (MSG_START, MSG_END, length)
    getFloatMeasurements();
    encodeAndSend(fDataOutput[0].asBytes, 4);
    encodeAndSend(fDataOutput[1].asBytes, 4);
    encodeAndSend(fDataOutput[2].asBytes, 4);
    resetMeasurementData();
    currentState = STATE_IDLE;
}


/** Handler of STATE_AUTOGAIN_ADCS */
void findOptimalADCGains(){
    /**Handler of STATE_AUTOGAIN_ADCS.

    Find the optimal gain for all ADCs.

    (1) measure several values. Remain in STATE_AUTOGAIN_ADCS until done.
    (2) for each ADC, pick a gain G such that the worst-case scenario
        measurement is at least ADC_RANGE_THRESHOLD/ADC_POSITIVE_SATURATION
        of the range of values that can be measured with G. The worst-case
        scenario value is calculated from the minimum and maximum values
        measured among the points.

    Reads
    -----
    numMeasurementsToDo, numMeasurementsDone, maximumPeak, minimumPeak,
    hardwareDetected, adc0Gain, adc1Gain

    Writes
    ------
    numMeasurementsDone, adc0Gain, adc1Gain

    Msg to PC
    ---------
    PC_OK : before returning to STATE_IDLE

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_AUTOGAIN_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If it takes more than 5s to acquire all the measurements
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_AUTOGAIN_ADCS (stays)
        While it has not yet finished measuring all the values required
    STATE_IDLE
        Successfully finished
    **/
    if (notInStateAndError(STATE_AUTOGAIN_ADCS))
        return;
    
    // Notice that we are not triggering the ADCs here, as we are
    // not particularly interested in a measuring from a specific
    // point in time. The first data will anyway be available
    // approximately 3/updateRate (i.e., 6 ms @ 500 Hz) after the
    // end of the self-calibration.
    measureADCsRipple();
    if (numMeasurementsDone < numMeasurementsToDo){
        if(didTimeOut()){
            // Do some cleanup:
            resetMeasurementData();
            numMeasurementsToDo = numMeasurementsToDoBackup;
            newMessage = false;

            // Place the ADCs back at gain zero, and
            // restore updateRate and calibration
            adc0Gain = 0;
            adc1Gain = 0;
            setAllADCgainsAndCalibration();
            if (currentState == STATE_ERROR)   // Some channel was not calibrated
                return;
        }
        return;
    }

    int16_t autogain_value0;
    int16_t autogain_value1;
    autogain_value0 = (max(abs(maximumPeak[0]), abs(minimumPeak[0]))
                       + (maximumPeak[0] - minimumPeak[0]));
    autogain_value1 = (max(abs(maximumPeak[1]), abs(minimumPeak[1]))
                       + (maximumPeak[1] - minimumPeak[1]));
    // TODO: probably something to check here: if either autogain_value is
    //       already in saturation with gain=0 something is wrong with
    //       the hardware or with the connections. Probably to check only
    //       if the specific ADC is present.
    // TODO: here we have to compute and store the ripple
    if(hardwareDetected.asInt & ADC_0_PRESENT){
        while(((autogain_value0 << (adc0Gain + 1)) < ADC_RANGE_THRESHOLD)
              && (adc0Gain < AD7705_MAX_GAIN)){
            adc0Gain++;
        }
    }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        while(((autogain_value1 << (adc1Gain + 1)) < ADC_RANGE_THRESHOLD)
              && (adc1Gain < AD7705_MAX_GAIN)){
            adc1Gain++;
      }
    }

    // Done. Clean up and update the ADCs with the new gain
    // found, also setting back the normal updateRate, and
    // loading the relevant calibration data.
    resetMeasurementData();
    numMeasurementsToDo = numMeasurementsToDoBackup;
    newMessage = false;
    setAllADCgainsAndCalibration();
    if (currentState == STATE_ERROR)   // Some channel was not calibrated
        return;

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}


/** Handler of STATE_ERROR */
void handleErrors(){
    /**Clean up after an error, and report it to the PC.
    
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
    if (notInStateAndError(STATE_ERROR))
        return;
    
    // First, report the error, so the PC knows
    // there may be some cleanup going on
    encodeAndSend(PC_ERROR);
    encodeAndSend(errorTraceback, 2);
    
    // Then clean up possible mess that caused the error
    switch(errorTraceback[1]){
        case ERROR_SERIAL_OVERFLOW:
            // Discard all characters in the serial buffer,
            // since messages are anyway likely corrupt.
            while(Serial.available()) Serial.read();
            break;
    }
    currentState = STATE_IDLE;
}


/** Handler of PC_RESET */
void reset(){
    /**
    Reset the Arduino, the ADCs and the DACs to the default
    settings, mimicking as much as possible the state right
    after boot-up or a hardware reset. The DAC output is
    also set to zero.

    Writes
    ------
    numBytesRead, msgLength, readingFromSerial, newMessage,
    waitingForDataFromPC, dacSettlingTime, hardwareDetected,
    adcUpdateRate, adc0Channel, adc1Channel, adc0Gain,
    adc1Gain, calibrationGain, adc0RipplePP, adc1RipplePP,
    adc0ShouldDecreaseGain, adc1ShouldDecreaseGain,
    numMeasurementsToDo, numMeasurementsToDoBackup,
    numMeasurementsDone, summedMeasurements, calibratedChannels
    
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
    dacSettlingTime = 100;

    hardwareDetected.asInt = 0;
    adcUpdateRate = AD7705_50HZ;
    adc0Channel = AD7705_CH0;
    adc1Channel = AD7705_CH0;
    adc0Gain = 0;
    adc1Gain = 0;
    calibrationGain = 0;
    adc0RipplePP = 0;
    adc1RipplePP = 0;
    adc0ShouldDecreaseGain = false;
    adc1ShouldDecreaseGain = false;

    numMeasurementsToDo = 1;
    numMeasurementsToDoBackup = 1;
    resetMeasurementData();

    // Rather than clearing the calibration data, we can just
    // reset the flags marking whether channels were calibrated
    // since we always check that they were calibrated before
    // loading calibration data
    for (byte iADC=0; iADC < N_MAX_ADCS_ON_PCB; iADC++){
        calibratedChannels[iADC][0] = false;
        calibratedChannels[iADC][1] = false;
    }

    // Explicitly reset the ADC I/O communication
    AD7705resetCommunication(CS_ADC_0);
    AD7705resetCommunication(CS_ADC_1);
    
    // Reset the update rate
    AD7705setClock(CS_ADC_0, adcUpdateRate);
    AD7705setClock(CS_ADC_1, adcUpdateRate);
    
    // As well as channel and gain
    AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);
    
    // Reset DAC and set the output voltage to zero
    AD5683reset(CS_DAC);
    AD5683setVoltage(CS_DAC, 0x0000);
}






/** ---------------- FUNCTIONS CALLED WITHIN STATE HANDLERS ---------------- **/

uint16_t getHardwarePresent() {
    /**
    Identify which hardware configuration we have on the board.

    The return value should be bitwise AND-ed with the *_PRESENT
    (for ADCs, LM35 and relay) and *_CLOSED (jumpers only) constants

    Hardware checks: ADCs, LM35, relay, jumpers

    Returns
    -------
    2-byte bit mask in a uint16_t
    */
    uint16_t result = 0;

    // (1) Check for ADCs
    delay(10);    // Make sure that input lines have settled (10 ms)

    // Trying to read a register (communication, but could
    // be anything else) will give 0xff if the line is held
    // high by the pull-ups on the MISO line or by the bus
    // isolator. This means the ADC is not there, or it is not
    // powered up (e.g., the external power supply is unplugged)
    byte adc0comm = AD7705readCommRegister(CS_ADC_0, AD7705_CH0);
    if (adc0comm != 0xff)
        result |= ADC_0_PRESENT;
    byte adc1comm = AD7705readCommRegister(CS_ADC_1, AD7705_CH0);
    if (adc1comm != 0xff)
        result |= ADC_1_PRESENT;

    // (2) Check for LM35 temperature sensor: the analog voltage
    // should be within the ADC range, and it should settle back
    // to a similar value after shortly connecting to a pull-up
    // resistor, then back to the LM35. In fact, since the LM35
    // cannot sink more than about 1 uA, connecting to a pull-up
    // resistor would drive the value high. Then, the value will
    // stay high if there is no LM35 there.

    // (2.1) Value before
    pinMode(LM35_PIN, INPUT);
    delay(10);    // Make sure the voltage has settled (10 ms)
    uint16_t sensorValue1 = analogReadMedian(LM35_PIN);

    // (2.2) Set pull-up resistor, and apply for 10 ms
    pinMode(LM35_PIN, INPUT_PULLUP);
    delay(10);

    // (2.3) Reset for normal measurements, and measure again
    pinMode(LM35_PIN, INPUT);
    delay(10);    // Make sure the voltage has settled (10 ms)
    uint16_t sensorValue2 = analogReadMedian(LM35_PIN);

    // (2.4) check the values
    if (sensorValue1 > 0 && sensorValue1 < LM35_MAX_ADU &&
        sensorValue2 < LM35_MAX_ADU && abs(sensorValue2 - sensorValue1) < 10)
        result |= LM35_PRESENT;

    // (3) Check for relay present: if the relay is mounted, there
    // should also be an external pull-up resistor soldered. This
    // gives about 0.12 V at the pin, i.e. 48 ADUs, i.e., in between
    // RELAY_MIN_ADU and RELAY_MAX_ADU
    uint16_t sensorValue = analogReadMedian(RELAY_PIN);
    // if (sensorValue > RELAY_MIN_ADU && sensorValue < RELAY_MIN_ADU)          // TODO: probably a bug here?? this goes well only if sensorValue == RELAY_MIN_ADU
    if (sensorValue > RELAY_MIN_ADU && sensorValue < RELAY_MAX_ADU)             // TODO: check correct: Replaced the second RELAY_MIN_ADU with RELAY_MAX_ADU
        result |= RELAY_PRESENT;
    else {
        // (4) Check jumper at JP3, which indicates that the user
        //     solidly selected the 2.5 V I0 range, and is not using
        //     a relay to switch ranges
        pinMode(JP_I0_PIN, INPUT_PULLUP);
        delay(1);
        if (digitalRead(JP_I0_PIN) == 0)
            result |= JP_I0_CLOSED;
        pinMode(JP_I0_PIN, INPUT);      // pull-up off, to reduce power consumption
    }

    // (5) Check jumper JP5 that should be closed by the
    // user if the 0--2.5 V AUX range as been chosen (at
    // JP6) instead of the 0--10 V AUX range
    pinMode(JP_AUX_PIN, INPUT_PULLUP);
    delay(1);
    if (digitalRead(JP_AUX_PIN) == 0)
        result |= JP_AUX_CLOSED;
    pinMode(JP_AUX_PIN, INPUT);
    return result;
}


bool isAcceptableADCUpdateRate(byte updateRate){                                // TODO: may later move to the ADC 'library'
    /**Return whether the given update rate in an acceptable value.*/
    return ((updateRate == AD7705_50HZ)
            or (updateRate == AD7705_60HZ)
            or (updateRate == AD7705_500HZ));
}


void selfCalibrateAllADCs(byte updateRate) {
    /**
    Simultaneous self-calibration for both ADCs
    (if present), using the current channel and gain.

    Parameters
    ----------
    updateRate : {AD7705_50HZ, AD7705_60HZ, AD7705_500HZ}
        The update rate for which calibration has to be performed
        Currently, no check is performed that the update rate
        is actually one of the acceptable values

    Reads
    -----
    adc0Channel, adc0Gain, adc1Channel, adc1Gain
    **/
    if (not isAcceptableADCUpdateRate(updateRate)){
        errorTraceback[0] = currentState;
        errorTraceback[1] = ERROR_RUNTIME;
        currentState = STATE_ERROR;
        return;
    }
    
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705selfCalibrate(CS_ADC_0, adc0Channel, adc0Gain, updateRate);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705selfCalibrate(CS_ADC_1, adc1Channel, adc1Gain, updateRate);
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705waitForCalibration(CS_ADC_0, adc0Channel);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705waitForCalibration(CS_ADC_1, adc1Channel);

    // AD7705 needs about 100--200 us to process the
    // self-calibration result. 1 ms delay is on the
    // safe side, and does not make much of a difference
    // on the overall procedure, which takes ~6/updateRate
    delay(1);
}


void setAllADCgainsAndCalibration() {
    /**
    Set the ADC gains according to the global adc0Gain,
    adc1Gain variables and restores the previously stored
    calibration values from the selfCalDataVsGain array.
    Also makes sure that the update rate of the ADCs is
    the correct one.

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel, adc0Gain,
    adc1Gain, selfCalDataVsGain, adcUpdateRate
    
    Goes to state
    -------------
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels was not calibrated before
    Same state
        Otherwise
    **/
    for (int iADC = 0; iADC < 2; iADC++) {
        // byte adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;       // BUG: ADC_0_PRESENT/ADC_1_PRESENT are uint16_t, not byte (probably would work with byte, but better have the right type to start with)!
        uint16_t adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;      // TODO: easier with externalADCs struct: externalADCs[i].present
        if (hardwareDetected.asInt & adcPresentMask) {
            byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;                    // TODO: easier with externalADCs[iADC].chipSelect
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: easier with externalADCs[iADC].channel
            byte gain = iADC==0 ? adc0Gain : adc1Gain;                          // TODO: easier with externalADCs[iADC].gain
            
            // Make sure that the channel has been calibrated before
            if (not calibratedChannels[iADC][channel]){
                errorTraceback[0] = currentState;
                errorTraceback[1] = ERROR_NEVER_CALIBRATED;
                currentState = STATE_ERROR;
                return;
            }
            
            int32_t cOffs = selfCalDataVsGain[gain][iADC][channel][0];
            int32_t cGain = selfCalDataVsGain[gain][iADC][channel][1];
            AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET, cOffs);
            AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN, cGain);
            
            // Reset the update rate. This is necessary after
            // finishing the auto-gain (that is done at a
            // different frequency), but does not hurt anyway
            AD7705setClock(chipSelect, adcUpdateRate);

            // Trigger to make sure the ADC is up to date with the new
            // calibration data, but we will not read the converted data
            AD7705setGainAndTrigger(chipSelect, channel, gain);
        }
    }
}


// Call this in other measurement functions that we may implement later on!
void decreaseADCGainsIfNeeded(){
    /**
    Decrease the gain of all ADCs if, during the last
    measurement, the value converted by the ADC has
    reached the upper part of the range.

    Writes
    ------
    adc0Gain, adc1Gain, adc0ShouldDecreaseGain, adc1ShouldDecreaseGain
    
    Goes to state
    -------------
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    Stays the same
        Otherwise
    **/
    if(adc0ShouldDecreaseGain){
        adc0Gain--;
        setAllADCgainsAndCalibration();    // Load the correct calibration
        if (currentState == STATE_ERROR)   // Channel was not calibrated
            return;
        adc0ShouldDecreaseGain = false;
    }
    if(adc1ShouldDecreaseGain){
        adc1Gain--;
        setAllADCgainsAndCalibration();    // Load the correct calibration
        if (currentState == STATE_ERROR)   // Channel was not calibrated
            return;
        adc1ShouldDecreaseGain = false;
    }
}


void makeAndSumMeasurements() {
    /**
    Get one measurement value from all the available ADCs, including the
    LM35, and add them up into summedMeasurements[]. The 'measurements'
    for devices not found are set to 0.

    Whenever a value is measured, it is checked against the saturation
    thresholds, possibly determining a gain switch the next time a
    measurement is performed.

    Notice that this does not trigger the ADCs, that should already be
    ready for measurement before this function is called.
    ADCs have to be triggered after changing the channel.

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel

    Writes
    ------
    summedMeasurements, numMeasurementsDone
    (indirectly) adc0Gain, adc1Gain, adc0ShouldDecreaseGain,
                 adc1ShouldDecreaseGain

    Msg to PC
    ---------
    None.

    Goes to state
    -------------
    STATE_ERROR : ERROR_ADC_SATURATED
        If one of the values read by the ADCs reaches a solid
        saturation, and the gain cannot be decreased further
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_MEASURE_ADCS (stays)
        Otherwise
    */
    // TODO: we should discuss what to do when a measurement saturates.
    //   Currently we are adding the measured value anyway. Perhaps we
    //   should shift the summedMeasurements by 1 and not increase the
    //   number of measurements? We could detect this by having
    //   checkMeasurementInADCRange return a boolean (true for in range,
    //   false for out of range), and add the measurement only if it is
    //   in range. In that case, we should take care of skipping completely
    //   this iteration for all the ADCs (i.e. not adding any of the
    //   measurements to summedMeasurements), and not increment the number
    //   of measurements at all (otherwise we screw the averaging).

    // TODO: we are not checking whether this value saturates (to
    //   LM35_MAX_ADU or an appropriate fraction of it, to be discussed),
    //   which may mean that there is a malfunction (e.g., high
    //   temperature). We may want to send an ERROR_TOO_HOT.

    int16_t measurement;
    if (hardwareDetected.asInt & LM35_PRESENT){
        measurement = analogReadMedian(LM35_PIN);
        summedMeasurements[2] += measurement;
        }   // this one first while we probably have to wait for the others
    if (hardwareDetected.asInt & ADC_0_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
        checkMeasurementInADCRange(adc0Gain, &adc0ShouldDecreaseGain,
                                   measurement);
        summedMeasurements[0] += measurement;
        }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        checkMeasurementInADCRange(adc1Gain, &adc1ShouldDecreaseGain,
                                   measurement);
        summedMeasurements[1] += measurement;
        }
    numMeasurementsDone++;
}


void checkMeasurementInADCRange(byte gain, bool* adcShouldDecreaseGain,
                                int16_t adcValue){    // The call to this function would be much easier having externalADCs, as then we would just pass the index of the ADC in the externalADCs array
    /**
    Check whether the (signed) value read is approaching
    the saturation for the current measurement range.

    Parameters
    ----------
    gain : byte
        The current gain of the ADC
    adcShouldDecreaseGain : *bool
        Pointer to the flag of the current ADC, which will be
        used to signal, in the next call to setVoltage(), that
        a gain reduction is necessary.
    adcValue : int16_t
        Signed version of the value measured by the ADC that
        has to be checked for saturation.
        (1) If adcValue is in the middle ~50% band of the
            measurement range (i.e., 2*ADC_RANGE_THRESHOLD),
            the value is still acceptable, but gain will soon
            need to be decreased. Set adcShouldDecreaseGain.
        (2) If adcValue is in solid saturation (full scale),                    // TODO: can't we avoid redoing all measurements? Perhaps it's enough to shift the correct summedMeasurements[i] by >> 1?
            the measurement itself is incorrect, and needs to
            be repeated from scratch after decreasing the gain.
        (3) If adcValue is in solid saturation but the gain
            cannot be decreased, go to STATE_ERROR,
            with ERROR_ADC_SATURATED code.

    Goes to state
    -------------
    STATE_ERROR : ERROR_ADC_SATURATED
        If adcValue is in solid saturation, and the gain cannot
        be decreased. Otherwise, the currentState is preserved.
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If the ADC channel used was not calibrated
    Stays unchanged
        Otherwise
    */
    // TODO: The ripple (measured at gain 0) should be >> by gain and
    //       added to abs(adcValue) for the next check only
    if(abs(adcValue) > ADC_RANGE_THRESHOLD
       && (gain > 0)
       && !*adcShouldDecreaseGain){
        // The measured value is above the "saturation" threshold,
        // but not yet at true saturation, which would make the
        // measured value completely wrong. Defer the decrease of
        // gain to someone else: currently this is done only when
        // setVoltage() gets called, i.e. at the next energy step.

        // NB: using the absolute value of the signed measurement
        //     means that we consider the value 'above threshold'
        //     when its unsigned version is outside a band twice
        //     as large as the one used in findOptimalADCGains()
        *adcShouldDecreaseGain = true;
    }

    if(((adcValue^=8000) == ADC_POSITIVE_SATURATION)
       || ((adcValue^=8000) == ADC_NEGATIVE_SATURATION)){
        if(gain > 0){
            gain--;
            setAllADCgainsAndCalibration();
            if (currentState == STATE_ERROR)   // Some channel was not calibrated
                return;
            resetMeasurementData();
            *adcShouldDecreaseGain = false;
        }
        // if(gain==0){       // TODO: I think this was a bug: if the previous 'if' decreased the gain to zero, this would trigger an error, even with a potentially acceptable value (measured with the lower gain).
        else{
            errorTraceback[0] = currentState;
            errorTraceback[1] = ERROR_ADC_SATURATED;
            currentState = STATE_ERROR;
        }
    }
}


void resetMeasurementData() {
    /** Reset summedMeasurements and numMeasurementsDone to 0 **/
    for (int i = 0; i < LENGTH(summedMeasurements); i++)
        summedMeasurements[i] = 0;
    numMeasurementsDone = 0;
}


void getFloatMeasurements() {
    /**
    Convert the global summedMeasurements to float in physical units
    (Volts, uAmps, degC), setting the value in the fDataOutput array.

    Notice that the I0 value (measured at the LEED electronics) will
    be returned in Volts, thus, the conversion to uAmps is delegated
    to the PC.

    Reads
    -----
    hardwareDetected, summedMeasurements, voltsPerBit, numMeasurementsToDo,
    adc0Channel, adc0Gain, adc1Channel, adc1Gain

    Writes
    ------
    fDataOutput
    **/
    if (hardwareDetected.asInt & ADC_0_PRESENT) {
        fDataOutput[0].asFloat = summedMeasurements[0] * voltsPerBit[adc0Gain] / numMeasurementsDone;
        // ADC#0, channel 0: I0 input (Volts).
        // The actual voltage at the input can be in either 0-2.5 V
        // (jumper closed) or 0-10 V (jumper open) ranges. This means
        // scaling the vale by ADC_0_CH0_SCALE_JO if the jumper is open,
        // and doing nothing if it's closed (ADC range is already 0-2.5V)
        if (adc0Channel == AD7705_CH0) {
            if ((hardwareDetected.asInt & JP_I0_CLOSED) == 0)  // 0-10 V with jumper open
                fDataOutput[0].asFloat *= ADC_0_CH0_SCALE_JO;
        }
        // ADC#0, channel 1: high voltage (Volts)
        else
            fDataOutput[0].asFloat *= ADC_0_CH1_SCALE;
    }

    if (hardwareDetected.asInt & ADC_1_PRESENT) {
        fDataOutput[1].asFloat = summedMeasurements[1] * voltsPerBit[adc1Gain] / numMeasurementsDone;
        // ADC#1, channel 0: I0 at biased sample (microAmps)
        if (adc1Channel == AD7705_CH0) {
            fDataOutput[1].asFloat *= ADC_1_CH0_SCALE;
        }
        else{
            // ADC#1 channel 1: AUX
            // The actual voltage at the input can be in either 0-2.5 V
            // (jumper closed) or 0-10 V (jumper open) ranges. This means
            // scaling the vale by JP_AUX_CLOSED if the jumper is open,
            // and doing nothing if it's closed (ADC range is already 0-2.5V) 
            if ((hardwareDetected.asInt & JP_AUX_CLOSED) == 0)  // 0-10 V with jumper open
                fDataOutput[1].asFloat *= ADC_1_CH1_SCALE_JO;
        }
    }

    if (hardwareDetected.asInt & LM35_PRESENT) // LM35: degrees C
        fDataOutput[2].asFloat = summedMeasurements[2] * ARDUINO_ADU_DEG_C / numMeasurementsDone;
}


void measureADCsRipple(){
    /**
    Measure a single value from any of the ADCs present, and
    update the maximum and minimum values measured so far.

    Reads
    -----
    hardwareDetected

    Writes
    ------
    maximumPeak, minimumPeak, numMeasurementsDone
    **/
    int16_t measurement;
    if (hardwareDetected.asInt & ADC_0_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
        if(measurement > maximumPeak[0]){
            maximumPeak[0] = measurement;
            }
        if(measurement < minimumPeak[0]){
            minimumPeak[0] = measurement;
            }
        }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        if(measurement > maximumPeak[1]){
            maximumPeak[1] = measurement;
            }
        if(measurement < minimumPeak[1]){
            minimumPeak[1] = measurement;
            }
        }
    numMeasurementsDone++;
}






/** -------------------------- ARDUINO UTILITIES --------------------------- **/

void setChipSelectHigh(byte ioPin) {
    /**
    Set a digital output used as a chip select signal from
    the default high-impedance state to high (=unselected),
    without a glitch to the low state.
    **/
    pinMode(ioPin, INPUT_PULLUP);
    digitalWrite(ioPin, HIGH);
    pinMode(ioPin, OUTPUT);
}






/** ---------------------------- OTHER UTILITIES --------------------------- **/

uint16_t analogReadMedian(byte pin) {
  /**
  Read the Arduino ADC for a given pin
  (A0...) three times and return the median
  */
  uint16_t a0 = analogRead(pin);
  uint16_t a1 = analogRead(pin);
  uint16_t a2 = analogRead(pin);
  return getMedian16(a0, a1, a2);
}

/* TODO: probably nicer to just have the getMedian, bigger and biggest
         functions just be overloaded for uint16_t and int32_t*/

/** Gets the median of three numbers */
uint16_t getMedian16(uint16_t a0, uint16_t a1, uint16_t a2) {
  uint16_t maximum = biggest16(a0, a1, a2);
  if (maximum == a0) return bigger16(a1, a2);
  if (maximum == a1) return bigger16(a0, a2);
  else return bigger16(a0, a1);
}

uint16_t bigger16(uint16_t a, uint16_t b) {
  return (a > b) ? a : b;
}

uint16_t biggest16(uint16_t a, uint16_t b, uint16_t c) {
  return bigger16(a, bigger16(b, c));
}

/** Gets the median of three numbers */
int32_t getMedian32(int32_t a0, int32_t a1, int32_t a2) {
  int32_t maximum = biggest32(a0, a1, a2);
  if (maximum == a0) return bigger32(a1, a2);
  if (maximum == a1) return bigger32(a0, a2);
  else return bigger32(a0, a1);
}

int32_t bigger32(int32_t a, int32_t b) {
  return (a > b) ? a : b;
}

int32_t biggest32(int32_t a, int32_t b, int32_t c) {
  return bigger32(a, bigger32(b, c));
}
