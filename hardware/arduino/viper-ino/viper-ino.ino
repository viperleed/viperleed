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
#define FIRMWARE_VERSION   0x0001    // MMmm, M=major, m=minor. Each of the four bytes can be 0--9 (MAX: 9999 == v99.99). CURENTLY: v0.1



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

    // Reset DAC
    AD5683reset(CS_DAC);

    // Explicitly reset the ADC I/O communication: just setting
    // the (inverted) chip select to high does not trigger a reset
    AD7705resetCommunication(CS_ADC_0);
    AD7705resetCommunication(CS_ADC_1);

    // Set DAC output to zero volts
    AD5683setVoltage(CS_DAC, 0x0000);

    #ifdef DEBUG
        Serial.print("hardware=0x");
        Serial.println(hardwareDetected.asInt, HEX);
    #endif
    #ifdef DEBUG
        pinMode(LED_BUILTIN, OUTPUT);
    #endif

    // Initialize values to their defaults
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
        case STATE_SETUP_ADC:
            setUpAllADCs();
            break;
        case STATE_SET_VOLTAGE:
            setVoltage();
            break;
        case STATE_TRIGGER_ADCS:
            if((millis() - initialTime) >= dacSettlingTime){
                triggerMeasurements();
            }
            break;
        case STATE_AUTOGAIN_ADCS:
            findOptimalADCGains();
            break;
        case STATE_MEASURE_ADCS:
            measureADCs();
            break;
        case STATE_ADC_VALUES_READY:
            sendMeasuredValues();
            break;
        case STATE_GET_HARDWARE:          // TODO: rename to STATE_GET_CONFIGURATION, as it also needs to return FIRMWARE_VERSION. Make the whole thing a single handler function like the others.
            hardwareDetected.asInt = getHardwarePresent();
            encodeAndSend(hardwareDetected.asBytes, 2);
            currentState = STATE_IDLE;
            break;
        case STATE_CALIBRATE_ADCS:
            calibrateADCsAtAllGains();
            break;
        case STATE_ERROR:                 // TODO: make a proper handler for this one!
            encodeAndSend(errorTraceback);
            currentState = STATE_IDLE;
            break;
    }
}

//================

void readFromSerial() {
    /**
    Receive one byte (i.e., character) from the serial line (i.e., PC) and store
    it into "temp_buffer[]". Also track whether all bytes that need to be
    received have been received by looking at whether the last byte read is
    MSG_END. When this happens, temp_buffer will be:
        [MSG_START, byte with length of message, message, MSG_END]

    Reads
    -----


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
        errorCode = ERROR_SERIAL_OVERFLOW;
        errorTraceback = currentState;
        currentState = STATE_ERROR;
        // TODO: read and discard whatever is in the buffer. This MUST be done
        //       in the error handler, where we will first send an error message
        //       and then clean up the buffer.
        return;
    }

    // TODO: We are currently reading only 1 'character' per sate loop          // ISSUE #11
    //       from the serial line. This means that it takes at least
    //       4 state-loop iterations to have the Arduino be responsive
    //       to a command. In the worst case, a single loop iteration
    //       currently takes ~360 ms (during STATE_CALIBRATE_ADCS).
    //       This means that it may take as long as ~1.4 seconds to
    //       Acknowledge a command and stop what's going on.
    //       Perhaps it would make more sense to read multiple characters
    //       from the serial line if there is anything available? In fact,
    //       if there are characters, it means that the PC sent some request
    //       that would anyway interrupt whatever we are doing after we read
    //       all the characters (one per state-loop iteration), so there is
    //       probably no point in waiting for whatever is going on to be over.

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
            errorCode = ERROR_MSG_TOO_LONG;
            errorTraceback = currentState;
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
    to a STATE_ERROR with the appropriate errorCode

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
        errorCode = ERROR_MSG_INCONSITENT;
        errorTraceback = currentState;
        currentState = STATE_ERROR;
        return false;
        }

    if (numDecodedBytes > 1){
        // Message is some data
        if (not waitingForDataFromPC){
            // But we're not expecting any
            errorCode = ERROR_MSG_UNKNOWN;
            errorTraceback = currentState;
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
        case PC_HARDWARE: break;
        case PC_INIT_ADC: break;
        case PC_MEASURE: break;
        case PC_RESET: break;
        case PC_SET_VOLTAGE: break;
        default:
            errorCode = ERROR_MSG_UNKNOWN;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            return false;
    }
    return true;
}

//============================
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
//============================

void updateState() {
/*
 * Use the message received from the python script to
 * update the state, provided that the Arduino is currently idle
 *
 * Parameters from globals:
 * ----------------
 * data_received[0]: byte array
 */
    if (not newMessage) return;

    // !!!!!!! need to remember to not always set newMessage = true             // TODO: What does this mean? when?
    switch(data_received[0]){
        case PC_HARDWARE:
            initialTime = millis();
            currentState = STATE_GET_HARDWARE;
            break;
        case PC_INIT_ADC:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SETUP_ADC;
            break;
        case PC_SET_VOLTAGE:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_VOLTAGE;
            break;
        case PC_AUTOGAIN:                                                       // TODO: may be nice to have a prepareForAutogain() function
            initialTime = millis();
            resetMeasurementData();
            adc0Gain = 0;
            adc1Gain = 0;
            selfCalibrateAllADCs(AD7705_500HZ);  // Takes ~13 ms
            numMeasurementsToDo = 25;                                           // TODO: back up numMeasurementsToDo before changing it!
            currentState = STATE_AUTOGAIN_ADCS;
            break;
        case PC_MEASURE:
            initialTime = millis();
            currentState = STATE_MEASURE_ADCS;
            break;
        case PC_CALIBRATION:
            // waitingForDataFromPC = true;  // TODO: will be the case after we rework this
            calibrationGain = 0;
            currentState = STATE_CALIBRATE_ADCS;
            break;
        case PC_RESET:
            reset();
            break;
    }
    newMessage = false;
}

//=========================

void findOptimalADCGains(){
    /**
    Find the optimal gain for all ADCs.

    (1) measure several values. Remain in AUTOGAIN_STATE until done.
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
    STATE_ERROR + ERROR_TIMEOUT : if it takes more than 5s to acquire
        all the measurements that are to be acquired
    STATE_AUTOGAIN_ADCS (stays) : while it has not yet finished measuring
        all the values required
    STATE_IDLE : successfully finished
    **/
    // Notice that we are not triggering the ADCs here, as we are
    // not particularly interested in a measuring from a specific
    // point in time. The first data will anyway be available
    // approximately 3/updateRate (i.e., 6 ms @ 500 Hz) after the
    // end of the self-calibration.
    measureADCsRipple();
    if (numMeasurementsDone < numMeasurementsToDo){
        if((millis() -  initialTime) > TIMEOUT){
            debugToPC("Timeout, Autogain failed, Arduino turns back to IDLE state");
            errorCode = ERROR_TIMEOUT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            resetMeasurementData();
            newMessage = false;
            // TODO: here we set gains to zero, and do the same as when going to IDLE at the end
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
    resetMeasurementData();
    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
    // TODO: here we have to do the same as in setUpAllADCs(), i.e. writing to updateRate the clock register, write the new gain to the ADCs, load the calibration, then return
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

//=========================
void setUpAllADCs(){
    /**
    Initialize the ADC from the parameters stored in the first 7 bytes of
    "data_received[]". The bytes have the following meaning:
    - [0] channel of ADC#0
    - [1] channel of ADC#1
    - [2] measurement frequency for ADCs #0 and #1
    - [3 + 4] number of measurements that will be averaged. [3] is the
              high byte, [4] the low byte. Data is interpreted as a uint16_t.
    - REMOVED [5] maximum gain value for ADC #0 and #1

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
    STATE_ERROR + ERROR_TIMEOUT : if more than 5s pass between the PC_INIT_ADC
        message and the receipt of data
    STATE_ERROR + ERROR_MSG_DATA_INVALID : if either channel or the update rate
        are not acceptable values
    STATE_SETUP_ADC (stays) : while waiting for data from the PC
    STATE_IDLE : successfully finished
    **/
    if(not newMessage){  // waiting for data from the PC
        if((millis() -  initialTime) > TIMEOUT){
            debugToPC("Timeout, no ini Values for ADC received!");              // TODO: probably comment out and have the Python take care of it
            errorCode = ERROR_TIMEOUT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            newMessage = false;
        }
        return;
    }
    
    waitingForDataFromPC = false;

    // TODO: we have to handle the bit mask here, so we also know if the
    // user wants us to measure the LM35 too.

    // Do some checking on the channels and updateRate values received:
    for (byte i = 0; i < 2; i++){
        byte channel = data_received[i];
        if (channel != AD7705_CH0 or channel != AD7705_CH1){
            errorCode = ERROR_MSG_DATA_INVALID;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            newMessage = false;
            return;
        }
    }

    adc0Channel = data_received[0];
    adc1Channel = data_received[1];

    adcUpdateRate = data_received[2];
    if (adcUpdateRate != AD7705_50HZ
            or adcUpdateRate != AD7705_60HZ
            or adcUpdateRate != AD7705_500HZ){
        errorCode = ERROR_MSG_DATA_INVALID;
        errorTraceback = currentState;
        currentState = STATE_ERROR;
        adcUpdateRate = AD7705_50HZ;
        newMessage = false;
        return;
    }

    numMeasurementsToDo = data_received[3] << 8 | data_received[4];

    //maximum_gain = data_received[5];                                          // TODO: remove from Python

    setAllADCgainsAndCalibration();                                             // TODO: this requires that the calibration has been properly done for the channels requested! We should check that this is in fact the case. This can perhaps be done in updateState().

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;//Back to command mode
    newMessage = false;
}
//=========================
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
    STATE_ERROR + ERROR_TIMEOUT : if more than 5s pass between the
        PC_SET_VOLTAGE message and the receipt of data
    STATE_SET_VOLTAGE (stays) : while waiting for data from the PC
    STATE_TRIGGER_ADCS : successfully finished
    **/
    if(not newMessage){   // Waiting for data from PC
        if((millis() -  initialTime) > TIMEOUT){
            debugToPC("Timeout, no ini Values for DAC received!");
            errorCode = ERROR_TIMEOUT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            newMessage = false;
        }
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
    if(adc0ShouldDecreaseGain){
        // During the last measurement, the value read by the ADC has reached the
        // upper part of the range. The gain needs to be decreased, and the ADC
        // requires recalibration
        adc0Gain--;
        setAllADCgainsAndCalibration();                                         // TODO: requires the calibration data to be available and up to date. This is currently not checked!
        adc0ShouldDecreaseGain = false;
    }
    if(adc1ShouldDecreaseGain){
        adc1Gain--;
        setAllADCgainsAndCalibration();
        adc0ShouldDecreaseGain = false;
    }

}

void measureADCs(){
    /**
    Acquire measurements from the ADCs.

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
    STATE_ERROR with ERROR_TIMEOUT : if it takes longer than 5s
        between the PC_MEASURE message and completing the number
        of measurements that need to be averaged
    STATE_ERROR with ERROR_ADC_SATURATED : if one of the inputs of
        the ADCs has reached saturation, and there is no room
        to decrease the gain.
    STATE_MEASURE_ADCS (stays) : until all the data values that
        need to be measured have been acquired
    STATE_ADC_VALUES_READY : successfully finished
    **/
    // TODO: since makeAndSumMeasurements() requires the ADCs to be already
    //       triggered fro measurements, the PC_MEASURE command cannot
    //       be handled until the ADCs have been triggered!
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

    if((millis() -  initialTime) > TIMEOUT){
        debugToPC("Timeout, ADC measurement Timeout!");
        errorCode = ERROR_TIMEOUT;
        errorTraceback = currentState;
        currentState = STATE_ERROR;
        newMessage = false;
    }
}
//=========================
void sendMeasuredValues(){
    /**
    Send measurements back to the PC.

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
    STATE_IDLE : always
    **/
    // TODO: Ideally, one would like to rather return a SINGLE MESSAGE
    //       containing at most 3*4 significant bytes (+ encoding), and
    //       including only the subset of measured values that were
    //       effectively requested. In this case, the worst-case scenario
    //       message length would be N_MAX_MEAS*4(bytes)*2(encoding)
    // TODO: we could check that currentState == STATE_ADC_VALUES_READY
    //       and throw an ERROR_RUNTIME otherwise.
    //       Possible drawback: if later we want to allow the PC to be
    //       able to force sending data even if not all the measurements
    //       have been completed, this would throw an unnecessary error
    getFloatMeasurements();
    encodeAndSend(fDataOutput[0].asBytes, 4); // !!!!!!!!!!!! need to implement way to only send demanded data
    encodeAndSend(fDataOutput[1].asBytes, 4);
    encodeAndSend(fDataOutput[2].asBytes, 4);
    resetMeasurementData();
    currentState = STATE_IDLE;
}

//=========================
void reset(){
    /**
    Reset the Arduino, the ADCs and the DACs to the default settings,
    i.e., the same as after boot-up or a hardware reset.

    The Arduino will be in STATE_IDLE at the end.
    **/
    currentState = STATE_IDLE;
    dacSettlingTime = 100;   // TODO: was 0, but should probably be the same as in the global initialization

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
    resetMeasurementData();

    // TODO: clear also calibration data

    AD7705resetCommunication(CS_ADC_0);
    AD7705resetCommunication(CS_ADC_1);
    AD5683reset(CS_DAC);
}

//=========================
void debugToPC(const char *debugmsg){  // !!!!!!!!!!!!!!! Since we're actually using this to send error messages back, I would rename it accordingly. Probably we should also define some error codes to accompany the message itself.
/*
 * Sends a debugmessage to the PC by masking the message
 * with a MSG_START +  DEBUGBYTE + the debug message +
 * MSG_END
 *
 * Parameters
 * ----------
 * debugmessage: const char array with variable field declaration
 */
  byte nb = PC_ERROR;
  Serial.write(MSG_START);
  Serial.write(nb);
  Serial.println(debugmsg);
  Serial.write(MSG_END);
}
//=========================
void debugToPC(byte debugbyte){
  byte nb = 0;
  Serial.write(MSG_START);
  Serial.write(nb);
  Serial.write(debugbyte);
  Serial.write(MSG_END);
}

void resetMeasurementData() {
    /** Reset summedMeasurement and numMeasurementsDone to 0 **/
    for (int i = 0; i < LENGTH(summedMeasurements); i++)
        summedMeasurements[i] = 0;
    numMeasurementsDone = 0;
}

void triggerMeasurements() {
    /**
    Triggers the ADCs to start the measurements now, for
    those ADCs that are available (i.e., AD7705_TRIGGER)

    Reads
    -----
    hardwareDetected, adc0Channel, adc0Gain, adc1Channel, adc1Gain

    Writes
    ------
    currentState

    Msg to PC
    ---------
    None.                                 // TODO: we actually have to return a PC_OK before going to STATE_MEASURE_ADCS

    Goes to state
    -------------
    STATE_IDLE : always                  // TODO: must always go to STATE_MEASURE_ADCS instead!
    **/
    // TODO: hardwareDetected needs to be a valid value before this
    // call makes sense at all. This is currently unchecked for, but
    // may be unnecessary (the function will do nothing if the PC did
    // send a PC_HARDWARE message before. What happens if hardwareDetected
    // is not up to date?
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);
    currentState = STATE_IDLE;
    // TODO: The ADCs are triggered for measurement, but NO MEASUREMENT IS
    //       SAVED! This means that the PC needs to request measurement
    //       explicitly. I don't really see then the use of the automatic
    //       switch from STATE_SET_VOLTAGE to STATE_TRIGGER_ADCS!!
}

void makeAndSumMeasurements() {
    /**
    Get one measurement value from all the available ADCs, including the
    LM35, and add them up into summedMeasurements[]. The 'measurements'
    for devices not found are set to 0.

    Whenever a value is measured, it is checked against the saturation
    thresholds, possibly determining a gain switch the next time a
    measurement is performed.                                                   // TODO: think if this is problematic when NOT running a ramp, but, say, measuring at constant voltage. Since the gain switch happens in setVoltage(), it would never get switched at all. This may be a problem when measuring transients.

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
    STATE_ERROR + ERROR_ADC_SATURATED : if one of the values read
        by the ADCs reaches a solid saturation, and the gain cannot
        be decreased to circumvent the problem.
    STATE_MEASURE_ADCS (stays) : unless a saturation error occurs
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
        checkMeasurementInADCRange(CS_ADC_0, adc0Channel, adc0Gain,
                                   &adc0ShouldDecreaseGain, measurement);
        summedMeasurements[0] += measurement;
        }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        checkMeasurementInADCRange(CS_ADC_1, adc1Channel, adc1Gain,
                                   &adc1ShouldDecreaseGain, measurement);
        summedMeasurements[1] += measurement;
        }
    numMeasurementsDone++;
}

void checkMeasurementInADCRange(byte chipSelectPin, byte channel, byte gain,
                                bool* adcShouldDecreaseGain, int16_t adcValue){    // The call to this function would be much easier having externalADCs, as then we would just pass the index of the ADC in the externalADCs array
    /**
    Check whether the (signed) value read is approaching
    the saturation for the current measurement range.

    Parameters
    ----------
    chipSelectPin : byte                                                        // TODO: unused
        The pin that can be used to communicate with the ADC
    channel : byte                                                              // TODO: unused
        The channel of the ADC that is currently being measured
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
    STATE_ERROR + ERROR_ADC_SATURATED : if adcValue is in solid
        saturation, and the gain cannot be decreased. Otherwise,
        the currentState is preserved.
    */
    // TODO: here we need to somehow use the adc0RipplePP and adc1RipplePP
    //       values, not just the bare adcValue measured! Alternatively, we
    //       could directly pass in the ripple-modified value during the call
    //         The ripple (measured at gain 0) should be >> by gain and
    //         added to abs(adcValue) for the next check only
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
            //AD7705selfCalibrate(chipSelectPin, channel, gain, adcUpdateRate);
            //AD7705waitForCalibration(chipSelectPin, channel);
            //delay(1);
            resetMeasurementData();
            *adcShouldDecreaseGain = false;
        }
        // if(gain==0){       // TODO: I think this was a bug: if the previous 'if' decreased the gain to zero, this would trigger an error, even with a potentially acceptable value (measured with the lower gain).
        else{
            debugToPC("ADC Overflow, maximum gain reached, no serious measurements possible...");
            errorCode = ERROR_ADC_SATURATED;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
        }
    }
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
        fDataOutput[0].asFloat = summedMeasurements[0] * voltsPerBit[adc0Gain] / numMeasurementsToDo;  // TODO: wouldn't it be better to use numMeasurementsDone?? They should be the same, but we may later want to allow the PC to send back data even if the averaging is not done, perhaps?
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
        fDataOutput[1].asFloat = summedMeasurements[1] * voltsPerBit[adc1Gain] / numMeasurementsToDo;
        // ADC#1, channel 0: I0 at biased sample (microAmps)
        if (adc1Channel == AD7705_CH0) {
            fDataOutput[1].asFloat *= ADC_1_CH0_SCALE;
        }
        else                            //ADC#1 channel 1: AUX
            if ((hardwareDetected.asInt & JP_AUX_CLOSED) == 0)  //0-10V with jumper open
                fDataOutput[1].asFloat *= ADC_1_CH1_SCALE_JO;
    }

    if (hardwareDetected.asInt & LM35_PRESENT) //LM35: degrees C
        fDataOutput[2].asFloat = summedMeasurements[2] * ARDUINO_ADU_DEG_C / numMeasurementsToDo;
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
        is actually one of the acceptable values                                // TODO: it may be a good idea to do this check

    Reads
    -----
    adc0Channel, adc0Gain, adc1Channel, adc1Gain
    **/
    // unsigned long startMillis = millis();                                    // TODO: remove, used only in (commented) DEBUG
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705selfCalibrate(CS_ADC_0, adc0Channel, adc0Gain, updateRate);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705selfCalibrate(CS_ADC_1, adc1Channel, adc1Gain, updateRate);
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705waitForCalibration(CS_ADC_0, adc0Channel);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705waitForCalibration(CS_ADC_1, adc1Channel);
    // unsigned long endMillis = millis();                                      // TODO: remove, used only in (commented) DEBUG
    //#ifdef DEBUG
    //    Serial.print("t_selfCal=");
    //    Serial.println(endMillis - startMillis);
    //#endif

    // AD7705 needs about 100--200 us to process the
    // self-calibration result. 1 ms delay is on the
    // safe side, and does not make much of a difference
    // on the overall procedure, which takes ~120 ms
    delay(1);
}

void storeAllSelfCalibrationResults(int32_t targetArray[2][2]) {                // TODO: use N_MAX_ADCS_ON_PCB for the dimension of the first direction.
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
    */
    for (int iADC = 0; iADC < 2; iADC++) {                                      // TODO: use N_MAX_ADCS_ON_PCB for the upper limit of the index
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

void setAllADCgainsAndCalibration() {
    /**
    Set the ADC gains according to the global adc0Gain,
    adc1Gain variables and restores the previously stored
    calibration values from the selfCalDataVsGain array.

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel, adc0Gain,
    adc1Gain, selfCalDataVsGain
    **/
    // TODO: requires the calibration data to be up to date (or at least
    //       that the self-calibration has run at least once)! This is
    //       currently unchecked for. Should probably raise some errors
    //       if they are not.
    for (int iADC = 0; iADC < 2; iADC++) {
        // byte adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;       // BUG: ADC_0_PRESENT/ADC_1_PRESENT are uint16_t, not byte (probably would work with byte, but better have the right type to start with)!
        uint16_t adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;      // TODO: easier with externalADCs struct: externalADCs[i].present
        if (hardwareDetected.asInt & adcPresentMask) {
            byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;                    // TODO: easier with externalADCs[iADC].chipSelect
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: easier with externalADCs[iADC].channel
            byte gain = iADC==0 ? adc0Gain : adc1Gain;                          // TODO: easier with externalADCs[iADC].gain
            int32_t cOffs = selfCalDataVsGain[gain][iADC][channel][0];
            int32_t cGain = selfCalDataVsGain[gain][iADC][channel][1];
            AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET, cOffs);
            AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN, cGain);
//          #ifdef DEBUG
//              Serial.print("gain=");
//              Serial.println(gain);
//              Serial.print(" cOffs=");
//              Serial.print(AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET));
//              Serial.print(" cGain=");
//              Serial.println(AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN));
//          #endif

            // Trigger to make sure the ADC is up to date with the new
            // calibration data, but we will not read the converted data
            AD7705setGainAndTrigger(chipSelect, channel, gain);
        }
    }
}

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
    // int result = 0;                                                          // TODO: probably was a bug: using int and returning uint16_t
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
    uint16_t sensorValue1 = analogReadMedian(LM35_PIN);                         // TODO: this and the following ones were int, but analogReadMedian returns uint16_t

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

    // (3) Check for relay present: if the relay is mounted, there              // TODO: I don't think we have anything implemented to actually drive the relay to switch gain. I can't recall if this is something that the user should just do externally with a switch, or if we wanted to allow switching this electronically
    // should also be an external pull-up resistor soldered. This
    // gives about 0.12 V at the pin, i.e. 48 ADUs, i.e. in between
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

void calibrateADCsAtAllGains(){
    /**
    For both ADCs, the channel currently selected is calibrated
    for all the possible gain values. The calibration results
    are stored, so that they can be fetched later for faster
    gain switching.

    Currently, each gain value is done over a separate state-loop.
    This means that it takes approximately 360 ms to finish a single
    run of this function (120 ms x 3 points for computing medians).
    This means, the serial line remains unresponsive for ~360 ms.
    TODO: We should consider whether we rather want to run one single
    median point per state-loop. This would leave the serial line
    inactive only for 120 ms.

    The calibration is done in parallel for both ADCs (if present).

    This function needs to be called again should the updateRate
    to be used for the measurements change. Since this takes quite
    long to complete (~3s for all the gains), it is advisable not
    to run it when finding the best ADC gain (don at 500Hz), since
    there the measurements are done at gain==0 and only for the
    channel to be measured.

    Writes
    ------
    calibrationGain, adc0Gain, adc1Gain

    Msg to PC
    ---------
    None.           // TODO: perhaps we would like to inform the PC once the whole calibration is done, since this is one of those things that takes time.


    Goes to state
    -------------
    STATE_CALIBRATE_ADCS (stays) : until all the gain values
        have been calibrated for the currently selected channels
    STATE_IDLE : after calibration is done
    */

    // TODO: I can see three issues with the current code:
    //   - it always does the calibration at AD7705_50HZ rather than at
    //     a specified updateRate. The problem is that the updateRate
    //     does not come as data from the PC until we go to STATE_SETUP_ADC.
    //     However, there we already need this calibration to be done.
    //     Solution: have the STATE_CALIBRATE_ADCS require a second
    //     communication from the PC (that may time out) in which we are
    //     told the updateRate to use. Then, we can get rid of the
    //     updateRate from the data required in STATE_SETUP_ADC.
    //   - A very similar problem exists for what concerns the channels:
    //     Here we need to know which channels we want to calibrate,
    //     but this information comes too late in STATE_SETUP_ADC,
    //     where we need the calibration already. Solution: have the
    //     STATE_CALIBRATE_ADCS require also the channels to be
    //     specified in the same communication message from the previous
    //     point. Rename this function fullyCalibrateOneADCChannel, and
    //     accept two parameters (adc0Channel and adc1Channel) that tell
    //     you which channels need to be calibrated. This way we can also
    //     allow the calibration to be done on both channels of both ADCs
    //     once we accept the bitmask (if three channels requested, it does
    //     not cost anything to do all four).
    //   - currently this function cannot time out, but it should! If the
    //     power plug for the ADCs is disconnected, we would wait forever
    //     for the calibration results (and the PC doesn't know how long
    //     it may take)!
    if (calibrationGain <= AD7705_MAX_GAIN) {
        adc0Gain = calibrationGain;
        adc1Gain = calibrationGain;
        // Calibrate 3 times each of the ADCs (in parallel) at
        // the current calibrationGain, and for the adc0Channel
        // and adc1Channel currently selected.                                  // TODO: in light of the comments above, one needs to set the channels passed as parameters at the beginning of this function
        // The three measurements are saved in selfCalDataForMedian,
        // which is later used to compute the 3-point median,
        // giving the final calibration result
        for (int i = 0; i < 3; i++) {
            selfCalibrateAllADCs(AD7705_50HZ);
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

        for (int iADC=0; iADC < 2; iADC++) {                                    // TODO: use N_MAX_ADCS_ON_PCB instead of 2
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
        // Calibration is over. Set gains to the lowest value and finish.
        adc0Gain = 0;
        adc1Gain = 0;
        currentState = STATE_IDLE;
    }
}


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




/* ---------- UTILITY FUNCTIONS ---------- */

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
