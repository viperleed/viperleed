/*
ViPErLEED - Firmware for Arduino hardware controller
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 16.04.2021
---------------------
*/

// ############################### HEADER START ################################

// Libraries
#include <Arduino.h>
#include <SPI.h>

// File with Arduino settings                                                   // TODO: make a header file
#include "viper-ino.h"

#define DEBUG              true      // Debug mode, writes to serial line, for use in serial monitor
#define FIRMWARE_VERSION   0x0001    // MMmm, M=major, m=minor. Each of the four can be 0--9 (MAX: 9999 == v99.99). CURENTLY: v0.1

// Useful macros, structs and unions
#define LENGTH(array)  (sizeof(array) / sizeof((array)[0]))
#ifndef NAN
    #define NAN   (1.0/0)   // Floating-point not-a-number
#endif

// Interpret the same value as either uint16_t or a 2-bytes array
union uint16OrBytes{
    uint16_t asInt;
    byte asBytes[2];
};

// Interpret the same value as either (32-bit) float or a 4-bytes array
union floatOrBytes{
    float asFloat;
    byte asBytes[4];
};




/** ------------------------- Communication with PC ------------------------ **/

// Constants for communication with the PC
#define MSG_START 254              // Beginning of a serial message
#define MSG_END 255                // End of a serial message
#define MSG_SPECIAL_BYTE 252       // Prevents clashing of data with (start, end, error) // TODO: CHANGE FOR PYTHON (Python has value 253)
#define MSG_MAX_LENGTH 16          // Max no. of bytes in an incoming serial message     // TODO: NEED TO INCREASE THIS TO HOST THE 3 FLOAT ADC VALUES
#ifndef SERIAL_BUFFER_SIZE
    #define SERIAL_BUFFER_SIZE 64  // Arduino limit for serial buffer. If buffer is full, new bytes are DISCARDED
#endif

// Acceptable messages for communication with the PC
#define PC_AUTOGAIN    8    // PC requested auto-gain for ADCs
#define PC_CALIBRATION 9    // PC requested self-calibration
#define PC_DEBUG       0                                                        // TODO: How is this different from PC_ERROR?
#define PC_ERROR     253    // An error occurred                                // TODO: Define in Python
#define PC_HARDWARE    3    // PC requested hardware configuration              // TODO: rename to make it clear that also firmware will be returned. Probably PC_GET_CONFIGURATION or similar
#define PC_INIT_ADC    4    // PC requested initialization of ADCs              // TODO: rename  // TODO: requires calibration data up to date as well as correct gain
#define PC_MEASURE     6    // PC requested to perform a measurement
#define PC_OK          5    // Acknowledge request from PC
#define PC_RESET      82    // PC requested a global reset (ASCII 'R')          // TODO: Change for python (Python has value 68)
#define PC_SET_VOLTAGE 7    // PC requested to set a certain energy

// Error codes
#define ERROR_NO_ERROR            0   // No error
#define ERROR_SERIAL_OVERFLOW     1   // Hardware overflow of Arduino serial
#define ERROR_MSG_TOO_LONG        2   // Too many characters in message from PC
#define ERROR_MSG_INCONSITENT     3   // Message received from PC is inconsistent. Probably corrupt.
#define ERROR_MSG_UNKNOWN         4   // Unknown request from PC
#define ERROR_MSG_DATA_INVALID    5   // Request from PC contains invalid information
#define ERROR_REQUEST_TOO_EARLY   6   // Command cannot be run before another one has been successfully completed
#define ERROR_TIMEOUT             7   // Timed out while waiting for something
#define ERROR_ADC_SATURATED       8   // One of the ADC values reached saturation, and gain can't be decreased further
#define ERROR_TOO_HOT             9   // The temperature read by the LM35 is too high
#define ERROR_RUNTIME           255   // Some function has been called from an inappropriate state. This is to flag possible bugs for future development.
byte errorCode;                       // Keeps track of which error occurred
byte errorTraceback;                  // Keeps track of the state that produced the error

// Variables used while communicating with the PC
byte numBytesRead = 0;                   // Counter for no. of bytes received. numBytesRead may be larger than the number of true data bytes due to encoding.
byte msgLength = 0;                      // Number of bytes to be expected in message, 2nd byte of received message
byte serialInputBuffer[MSG_MAX_LENGTH];  // Contains all the raw (i.e., still encoded) bytes read from the serial line
boolean readingFromSerial = false;       // True while Arduino is receiving a message
boolean newMessage = false;              // True when a complete, acceptable message has been read

/** TODO
    Judging from the code, and if I didn't misinterpret anything, I think
    we can have only one buffer for the input/output to the PC, so we can
    replace "data_received" and "data_send" with a single "message" array
    Alternatively, we can have them be a decodedCommand[] and encodedReply[]
*/
byte data_received[MSG_MAX_LENGTH];      // Real message
byte data_send[MSG_MAX_LENGTH];




/** ------------------------- Finite state machine ------------------------- **/

#define STATE_IDLE                 0  // Wait for requests from PC
#define STATE_SETUP_ADC            1  // Pick correct ADC channels, update frequency, and no. of measurement points  // TODO: comment says it also triggers measurements!  // TODO: Maybe rename to STATE_PREPARE_ADCS_FOR_MEASUREMENT, STATE_PREPARE_FOR_MEASUREMENT, or STATE_SETUP_ADCS?
#define STATE_SET_VOLTAGE          2  // Set a voltage with the DAC
#define STATE_TRIGGER_ADCS         3  // Start a measurement right now
#define STATE_ADC_MEASURE          4  // ADC measurements in progress
#define STATE_ADC_VALUES_READY     5  // ADC measurements done
#define STATE_AUTOGAIN             6  // Find optimal gain for both ADCs
#define STATE_GET_HARDWARE         7  // Find current hardware configuration
#define STATE_INITIAL_CALIBRATION  8  // Figure out correct offset and calibration factors for ADCs at all gains.
#define STATE_ERROR                9  // An error occurred
uint16_t currentState = STATE_IDLE;




/** ---------------------- Hardware-specific settings ---------------------- **/
// Measurement devices
#define N_MAX_MEAS         3  // Number of measurement devices that we can have: ADC#0, ADC#1, LM35 (if present)
#define N_MAX_ADCS_ON_PCB  2  // Maximum number of ADCs on the PCB
byte numADCsOnBoard = N_MAX_ADCS_ON_PCB;  // Real number of ADCs present      TODO: have getHardwarePresent set this value // TODO: I forgot why I actually wanted to use this. Probably easier to just store as many as N_MAX_ADCS_ON_PCB

// I/O pins on Arduino board
#define CS_DAC             4    // (Inverted) chip select of DAC
#define CS_ADC_0           5    // (Inverted) chip select of ADC #0
#define CS_ADC_1           6    // (Inverted) chip select of ADC #1 (optional)
#define LM35_PIN          A1    // LM35 temperature sensor (optional)
#define RELAY_PIN         A4    // Relay is connected here if present
#define JP_I0_PIN         A4    // Jumper for manual I0 2.5V range, same as relay pin
#define JP_AUX_PIN        A5    // Jumper for manual AUX 2.5 V range

// Arduino internal reference, ADC maximum, and settings for relay an LM35
#define VREF_INTERNAL   2.56    // Arduino micro internal ADC reference is 2.56 V
#define ARDUINO_ADC_MAX  0x3ff  // 10-bit internal ADC
#define RELAY_MIN_ADU     24    // Relay has about 0.12 V with pull-up, ADC signal is larger than this
#define RELAY_MAX_ADU     96    // ADC signal of relay with pull-up is less than this
#define LM35_MAX_ADU     320    // LM35 signal should be less than 0.8 V (80 degC)
#define ARDUINO_ADU_VOLTS (VREF_INTERNAL/ARDUINO_ADC_MAX)      // Internal ADC: volts per ADU
#define ARDUINO_ADU_DEG_C (100*VREF_INTERNAL/ARDUINO_ADC_MAX)  // LM35 degC per ADU

#define R_SOURCE        1300.0   // Source resistance of our circuit at ADC inputs, 1.3 kOhm. Used in voltsPerBit[] below.

// Scaling of ADC input ranges due to voltage dividers or preamplification on the board
#define ADC_0_CH0_SCALE_JO      0.004          // ADC#0 channel 0: with JP at I0 open or relay off, in volts
#define ADC_0_CH1_SCALE    (16000./39.*0.001)  // ADC#0 channel 1: High-voltage divider
#define ADC_1_CH0_SCALE        (10/2500.0)     // ADC#1 channel 0: uAmps
#define ADC_1_CH1_SCALE_JO      0.004          // ADC#1 channel 1: with JP5 open (voltage divider), in volts

// Which hardware and closed jumpers do we have? Bits of present/closed jumpers are 1
#define ADC_0_PRESENT    0x01 // Bit is set if ADC #0 was detected
#define ADC_1_PRESENT    0x02 // Bit is set if ADC #1 was detected
#define LM35_PRESENT     0x04 // Bit is set if LM35 temperature sensor was detected
#define RELAY_PRESENT    0x08 // Bit is set if the relay for I0 input 2.5V/10V is present
#define JP_I0_CLOSED     0x10 // Bit is set if JP3 7-8 is closed or relay on (to indicate 2.5V I0 range)
#define JP_AUX_CLOSED    0x20 // Bit is set if JP5 is closed (to indicate 2.5V AUX range rather than 10V range)

// ADC saturation thresholds
#define ADC_POSITIVE_SATURATION 0xffff  // Max ADC output in a range
#define ADC_NEGATIVE_SATURATION 0x0000  // Min ADC output in a range
#define ADC_RANGE_THRESHOLD     0x3fff  // If output is larger than this (~25% of range), we need a new gain next time




/** --------------------- Analog-to-digital converters ---------------------- */
/** NOTE:
    The contents of this section are very specific to the AD7705
    ADCs. Perhaps it would be nice to create a generic ADC library
    (or ADC class?) that can handle this stuff. The code would then
    look much cleaner.
*/

// Definitions for AD7705 analog-to-digital converter
#define AD7705_SPIMODE     3   // Data accepted at rising edge of SCLK
#define AD7705_MAX_GAIN    7   // Gains are 0 (x1) to 7 (x128)

//AD7705 communication register
#define AD7705_READ_REG   0x08  // Set bit for reading a register
#define AD7705_CH0        0x00  // Set bit for addressing channel 0
#define AD7705_CH1        0x01  // Set bit for addressing channel 1
#define AD7705_REG_COMM   0x00  // Bits to set for accessing communication register
#define AD7705_REG_SETUP  0x10  // Bits to set for accessing setup register
#define AD7705_REG_CLOCK  0x20  // Bits to set for accessing clock register
#define AD7705_REG_DATA   0x30  // Bits to set for accessing data register (16 bit)
#define AD7705_REG_OFFSET 0x60  // Bits to set for accessing calibration offset register (24 bit)
#define AD7705_REG_GAIN   0x70  // Bits to set for accessing calibration gain register (24 bit)
#define AD7705_DRDY       0x80  // Bit mask for Data Ready bit in communication register

//AD7705 setup register
#define AD7705_TRIGGER    0x01  // Bit to set for triggering ADC (called FSYNC in data sheet)
#define AD7705_BUFMODE    0x02  // Bit to set for buffered mode (not used by us)
#define AD7705_UNIPOLAR   0x04  // Bit to set for unipolar mode (not used by us)
#define AD7705_SELFCAL    0x40  // Bit to set for self-calibration

//AD7705 clock register
#define AD7705_CLK        0x0C  // Bits to set for 4.9152 MHz crystal: CLK=1, CLKDIV=1
#define AD7705_50HZ       0x04  // Bits to set for 50 Hz update rate and suppression
#define AD7705_60HZ       0x05  // Bits to set for 60 Hz update rate and suppression
#define AD7705_500HZ      0x07  // Bits to set for 500 Hz update rate

// Correction for the finite source impedance and ADC input impedance
#define R_IN_0              4.1e6       // input resistance of AD7705 in gain0 = x1, measured typical value
#define REF_OVER_RANGE  (2500.0/32768)  // millivolts (uA@1kOhm) per bit at gain0, bipolar, input 0ohm

const float voltsPerBit[] = {
    REF_OVER_RANGE/R_IN_0*(R_IN_0 + R_SOURCE),         // gain0 = x1
    REF_OVER_RANGE/R_IN_0*2*(R_IN_0/2 + R_SOURCE)/2,   // gain1 = x2 has half R_IN_0
    REF_OVER_RANGE/R_IN_0*4*(R_IN_0/4 + R_SOURCE)/4,   // gain2 = x4 has 1/4 R_IN_0
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/8,   // gain3 = x8 and up: 1/8 R_IN_0
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/16,  // gain4 = x16
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/32,  // gain5 = x32
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/64,  // gain6 = x64
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/128, // gain7 = x128
};




/** ---------------------- Digital-to-analog converter ---------------------- */
/**
TODO: This stuff will go into the DAC library
*/

// Definitions for AD5683 digital-to-analog converter
#define AD5683_SPIMODE   1    // Data accepted at falling edge of SCLK
#define AD5683_SET_DAC   0x30 // Bits of highest (first) byte for writing & setting DAC
                              // Lower 4 bits must be highest 4 bits of data
#define AD5683_RESET_MSB 0x40 // Bits of highest (first) byte for reset (internal reference)




/** ----------------------- ADC and DAC communication ----------------------- */
// SPI communication settings
SPISettings AD5683_SPI_SETTING(2500000, MSBFIRST, AD5683_SPIMODE);
SPISettings AD7705_SPI_SETTING(2500000, MSBFIRST, AD7705_SPIMODE);




/** -------------------- Globals for firmware functions -------------------- **/
// Timers (defined in milliseconds)
#define TIMEOUT 5000                  // Max 5 seconds to do stuff
unsigned long initialTime;            // System time when switching to a new state
uint16_t      dacSettlingTime = 100;  // The time interval for the DAC output to be stable
                                      //   This is just a default value. The actual one comes
                                      //   from the PC with a PC_SET_VOLTAGE command, and is
                                      //   read in setVoltage()

/* union integerOrBytes{
  uint16_t asInt;
  byte asBytes[2];
}             hardwareDetected;       // Bits set indicate this hardware is present/jumper closed */
uint16OrBytes hardwareDetected;       // Bits set indicate this hardware is present/jumper closed

// ADCs: measurement frequency, channel, gain
byte     adcUpdateRate;        // Update rate for both ADCs (will be set for line frequency)
byte     adc0Channel;          // Channel of ADC#0
byte     adc1Channel;          // Channel of ADC#1
byte     adc0Gain = 0;         // Gain of ADC#0, 0...7
byte     adc1Gain = 0;         // Gain of ADC#0, 0...7
bool     adc0ShouldDecreaseGain = false;  // Whether ADC#0 should increase its gain in the next cycle
bool     adc1ShouldDecreaseGain = false;  // Whether ADC#0 should increase its gain in the next cycle

// ADCs: quantities needed for self-calibration
byte     calibrationGain = 0;           // Gain for which a calibration is currently being performed (in parallel for both ADCs)
int32_t  selfCalDataForMedian[3][2][2]; // Values from which medians are calculated: three for median, two ADCs, last index is offset(0) & gain(1) // TODO: this is currently used only locally inside initialCalibration(). If we keep the code as is, this should move in there. If we decide to do one median per state-loop it can stay out here.
int32_t  selfCalDataVsGain[AD7705_MAX_GAIN + 1][2][2][2]; // For each gain, two ADCs, two channels each, and last index is offset(0)&gain(1)

// ADCs: variables for measuring and storing the line frequency ripple
int16_t  maximumPeak[N_MAX_ADCS_ON_PCB];       // Maximum of measurement, one for each ADC
int16_t  minimumPeak[N_MAX_ADCS_ON_PCB];       // Minimum of measurement, one for each ADC
uint16_t adc0RipplePP = 0;     // Ripple (peak-peak) measured at gain=0 for ADC#0 during auto-gain
uint16_t adc1RipplePP = 0;     // Ripple (peak-peak) measured at gain=0 for ADC#1 during auto-gain

/* // ADC container                                                             TODO: use it to simplify code below
struct analogToDigitalConverter {
    byte     chipSelect;    // to be initialized!
    uint16_t present;       // to be initialized!
    byte     channel = AD7705_CH0;
    byte     gain = 0;
    bool     shouldDecreaseGain = false;
    int16_t  rippleMinimum = 0;
    int16_t  rippleMaximum = 0;
    uint16_t ripplePP = 0;                                   // Ripple (peak-peak) measured at gain=0 during auto-gain
    int32_t  calibrationDataForMedian[3][2];                 // Three values for median, last index is offset(0) & gain(1) // TODO: perhaps this is not needed and can stay as it currently is, given the way it is currently implemented.
    int32_t  calibrationForGain[AD7705_MAX_GAIN + 1][2][2];  // For each gain, two channels each, and last index is offset(0) & gain(1)
} externalADCs[N_MAX_ADCS_ON_PCB]; */

// Measurements
uint16_t numMeasurementsToDo = 1;         // No. of ADC measurements to do before measurement is considered over
uint16_t numMeasurementsDone;             // Counter for the number of ADC measurements done so far
int32_t  summedMeasurements[N_MAX_MEAS];  // Measurements of ADC#0, ADC#1 and LM35 are summed up here

/*
union floatOrBytes{                       // Measured ADC voltages, as floats and
  float asFloat;                          // 4-byte array, useful for sending back
  byte asBytes[4];                        // to PC the values measured by the ADCs
} fDataOutput[N_MAX_MEAS]; */
floatOrBytes fDataOutput[N_MAX_MEAS];     // Measurements as voltages  // TODO: rename measuredVoltages[]

// ################################ HEADER END #################################




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
        case STATE_AUTOGAIN:
            findOptimalADCGains();
            break;
        case STATE_ADC_MEASURE:
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
        case STATE_INITIAL_CALIBRATION:   // TODO: rename
            initialCalibration();
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
    MSG_END. When this happens, temp_buffer will be:              !!!!!!!!!!!!!!!! Probably it would also be good to check here if the number of bytes read is the one expected?
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

    if(Serial.available() == SERIAL_BUFFER_SIZE){                               // TODO: Was '>=', but cannot be '>'
        // Serial.println("Serial overflow.");                                  // TODO: This is probably unnecessary if we send back also the errorCode
        errorCode = ERROR_SERIAL_OVERFLOW;
        errorTraceback = currentState;
        currentState = STATE_ERROR;
        // TODO: we're not doing anything to solve the overflow issue!!
        return;
    }

    // TODO: perhaps it would make more sense to read multiple characters
    //       from the serial line if there is anything available? In fact,
    //       if there are characters, it means that the PC sent some request
    //       that would anyway interrupt whatever we are doing after we read
    //       all the characters (one per state-loop iteration), so there is
    //       probably no point in waiting for whatever is going on to be over.
    
    // TODO: make sure we do not go into a segmentation fault by writing
    //       to places we shouldn't! In fact, our serialInputBuffer has
    //       a limited size of MSG_MAX_LENGTH bytes. Here we're not
    //       making in any way sure that we read at most MSG_MAX_LENGTH bytes.
    //       If we are going to, we should throw an ERROR_MSG_TOO_LONG, 
    //       and we should anyway read (and throw away) all the characters
    //       till the next MSG_END to prevent leaving a half-message hanging.
    //       This is a little tricky to do, however (it may be that part
    //       of the message has still to be delivered to the serial).
    //       I think we should anyway check the message received for
    //       consistency [e.g., a 'command' message should be of length
    //       exactly == 1(length)+1(the command); we do expect to receive
    //       data only when we are in certain states; and similar). Then
    //       we can throw an ERROR_MSG_UNKNOWN or ERROR_MSG_INCONSITENT
    //       to signal the issue. However, this means that an ERROR_MSG_TOO_LONG
    //       can be often followed by another error message. Another option
    //       is to not always immediately return to STATE_IDLE from
    //       STATE_ERROR: stay in STATE_ERROR if errorCode is ERROR_MSG_TOO_LONG
    //       until we receive a MSG_END character that marks the end of the
    //       too-long message to be discarded.
    byte byteRead = Serial.read();
    
    // New message
    if (byteRead == MSG_START) {
        numBytesRead = 0;
        readingFromSerial = true;
    }

    // Accumulate characters
    if(readingFromSerial) {
        serialInputBuffer[numBytesRead] = byteRead;
        numBytesRead++;
    }

    // A full message has been read
    if (byteRead == MSG_END) {
        readingFromSerial = false;
        msgLength = serialInputBuffer[1];
        if (msgLength != numBytesRead){                                         // TODO: check that this is correct. Depends on whether msg_length includes the MSG_END byte or not!
            errorCode = ERROR_MSG_INCONSITENT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            return;
        }
        newMessage = true;  // !!!!!! make sure to not always set to true       // TODO: What does this mean? When would you not set it to true?
        decodeMessage();
    }
}

void decodeMessage(){
/*
    Move the interesting bytes of serialInputBuffer[] into data_received[].
     *
    In practice, only the actual characters are kept. This means:
    (1) Skipping:
           MSG_START == serialInputBuffer[0]
           total no. of bytes in message == serialInputBuffer[1]
           MSG_END == serialInputBuffer[last]
    (2) Decoding bytes with value MSG_SPECIAL_BYTE. In this case, the
        actual character is MSG_SPECIAL_BYTE + the next character.

    Returns into globals
    --------------------
    data_received[] : byte array
*/
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
            initialTime = millis();
            currentState = STATE_SETUP_ADC;
            break;
        case PC_SET_VOLTAGE:
            initialTime = millis();
            currentState = STATE_SET_VOLTAGE;
            break;
        case PC_AUTOGAIN:                                                       // TODO: may be nice to have a prepareForAutogain() function
            initialTime = millis();
            resetMeasurementData();
            adc0Gain = 0;
            adc1Gain = 0;
            adcUpdateRate = AD7705_500HZ;
            selfCalibrateAllADCs(AD7705_500HZ);  // Takes ~13 ms
            numMeasurementsToDo = 20; // !!!!!!!!! May need to be changed later on  // TODO: If I recall correctly, we said 25 measurement points here, not 20!
            currentState = STATE_AUTOGAIN;
            break;
        case PC_MEASURE:
            initialTime = millis();
            currentState = STATE_ADC_MEASURE;
            break;
        case PC_CALIBRATION:
            calibrationGain = 0;
            currentState = STATE_INITIAL_CALIBRATION;
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
    STATE_AUTOGAIN (stays) : while it has not yet finished measuring
        all the values required
    STATE_IDLE : successfully finished
    **/
    // Notice that we are not triggering the ADCs here, as we are
    // not particularly interested in a measuring from a specific
    // point in time. The first data will anyway be available
    // approximately 3/updateRate (i.e., 6 ms @ 500 Hz) after the
    // end of the self-calibration.
    measureADCsAndPeaks();                                                      // TODO: we may be waiting forever if the ADCs power is disconnected during this! What happens in AD7705waitAndReadData in this case?
    if (numMeasurementsDone < numMeasurementsToDo){
        if((millis() -  initialTime) > TIMEOUT){
            debugToPC("Timeout, Autogain failed, Arduino turns back to IDLE state");
            errorCode = ERROR_TIMEOUT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            resetMeasurementData();
            newMessage = false;
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
}

void measureADCsAndPeaks(){
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
    if(not newMessage){
        if((millis() -  initialTime) > TIMEOUT){
            debugToPC("Timeout, no ini Values for ADC received!");              // TODO: probably comment out and have the Python take care of it
            errorCode = ERROR_TIMEOUT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            newMessage = false;
        }
        return;
    }

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

    //maximum_gain = data_received[5]; //!!!!!!!!!! remove from Python

    //selfCalibrateAllADCs(adcUpdateRate);
    setAllADCgainsAndCalibration();                                             //  TODO: this requires that the calibration has been properly done for the channels requested! We should check that this is in fact the case. This can perhaps be done in updateState().

    //delay(1);
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
    None. Probably intended not to cram the serial line when running
    a voltage ramp.

    Goes to state
    -------------
    STATE_ERROR + ERROR_TIMEOUT : if more than 5s pass between the
        PC_SET_VOLTAGE message and the receipt of data
    STATE_SET_VOLTAGE (stays) : while waiting for data from the PC
    STATE_TRIGGER_ADCS : successfully finished
    **/
    if(not newMessage){
        if((millis() -  initialTime) > TIMEOUT){
            debugToPC("Timeout, no ini Values for DAC received!");
            errorCode = ERROR_TIMEOUT;
            errorTraceback = currentState;
            currentState = STATE_ERROR;
            newMessage = false;
        }
        return;
    }

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
    STATE_ERROR + ERROR_TIMEOUT : if it takes longer than 5s
        between the PC_MEASURE message and completing the number
        of measurements that need to be averaged
    STATE_ERROR + ERROR_ADC_SATURATED : if one of the inputs of
        the ADCs has reached saturation, and there is no room
        to decrease the gain.
    STATE_ADC_MEASURE (stays) : until all the data values that
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
  byte nb = PC_DEBUG;
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
    None.

    Goes to state
    -------------
    STATE_IDLE : always
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
    STATE_ADC_MEASURE (stays) : unless a saturation error occurs
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

    Both ADCs are also triggered to start converting
    the analog values at adc0Channel and adc1Channel.

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

void initialCalibration(){      // TODO: rename
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
    calibrationGain, adc0Gain, adc1Gain, 

    Msg to PC
    ---------
    None.           // TODO: perhaps we would like to inform the PC once the whole calibration is done, since this is one of those things that takes time.
    

    Goes to state
    -------------
    STATE_INITIAL_CALIBRATION (stays) : until all the gain values
        have been calibrated for the currently selected channels
    STATE_IDLE : after calibration is done
    */

    // TODO: I can see three issues with the current code:
    //   - it always does the calibration at AD7705_50HZ rather than at
    //     a specified updateRate. The problem is that the updateRate
    //     does not come as data from the PC until we go to STATE_SETUP_ADC.
    //     However, there we already need this calibration to be done.
    //     Solution: have the STATE_INITIAL_CALIBRATION require a second
    //     communication from the PC (that may time out) in which we are
    //     told the updateRate to use. Then, we can get rid of the
    //     updateRate from the data required in STATE_SETUP_ADC.
    //   - A very similar problem exists for what concerns the channels:
    //     Here we need to know which channels we want to calibrate,
    //     but this information comes too late in STATE_SETUP_ADC,
    //     where we need the calibration already. Solution: have the
    //     STATE_INITIAL_CALIBRATION require also the channels to be
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


/* ---------- AD7705 ADC FUNCTIONS ---------- */

/** Starts an I/O operation for the AD7705 with given chip select pin */
void AD7705startIO(byte chipSelectPin) {
    SPI.beginTransaction(AD7705_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
}

/** Finishes an I/O operation for the AD7705 with given chip select pin */
void AD7705endIO(byte chipSelectPin) {
    digitalWrite(chipSelectPin, HIGH);
}

/** Resets I/O of the AD7705 by writing 32 high bits */
void AD7705resetCommunication(byte chipSelectPin) {
    AD7705startIO(chipSelectPin);
    SPI.transfer16(0xffff);
    SPI.transfer16(0xffff);
    AD7705endIO(chipSelectPin);
}

/** Writes the clock register of the AD7705 with the given chip select pin
 *  The update rate can be AD7705_50HZ, AD7705_60HZ, or AD7705_500HZ
 *  Not used because setting the update rate requires self-calibration */
 /*void AD7705setClock(byte chipSelectPin, byte updateRate) {
   AD7705startIO(chipSelectPin);
   SPI.transfer(AD7705_REG_CLOCK);
   SPI.transfer(AD7705_CLK | updateRate);
   AD7705endIO(chipSelectPin);
 }*/

void AD7705setGainAndTrigger(byte chipSelectPin, byte channel, byte gain) {
    /**Set the gain and trigger an AD7705 ADC.

    This essentially sets a well-defined point in time after
    which conversions of the selected input channel will be
    available, after amplification with the selected gain.
    
    Notice that this function does not initiate the reading
    of any of the conversion results from the Arduino. Thus,
    every 1/updateRate sec after this, new values will be
    available in the data register, and they will be thrown
    away if they are not read.
    
    The function also sets bipolar unbuffered operation mode,
    as well as 'normal operation' (i.e., not a calibration)
    that is necessary for our inputs. Any old data that may
    still be present in the data register is also discarded.

    Parameters
    ----------
    chipSelectPin : byte
        The Arduino pin that corresponds to the ADC that has to
        be addressed
    channel : byte, {0, 1}
        Which of the input channels will be used for conversion
    gain : byte, {0...7}
        Which gain value should be used. The input is amplified
        by a factor 2**gain.
    **/
    AD7705startIO(chipSelectPin);
    //SPI.transfer16(0xffff); //reset should not be necessary
    //SPI.transfer16(0xffff);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_TRIGGER | gain << 3);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(gain << 3);

    // Read the data register to discard any old data
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0xffff);
    AD7705endIO(chipSelectPin);
}

void AD7705selfCalibrate(byte chipSelectPin, byte channel,
                         byte gain, byte updateRate) {
    /**Self-calibrate an AD7705 ADC.
    
    Thereafter, one has to wait until self-calibration is
    done, by calling AD7705waitForCalibration, or by waiting
    for the first data (AD7705waitAndReadData).
    
    Triggering (STATE_TRIGGER_ADCS) during self-calibration
    should never be performed.
    
    Parameters
    ----------
    chipSelectPin : byte
        The Arduino pin corresponding to the ADC
        that needs to be calibrated.
    channel : byte, {0, 1}
        Which of the input channels needs to be
        calibrated.  Only this channel will be
        calibrated. The calibration of the other
        channel remains untouched.
    gain : byte, {0...7}
        Which gain should be used for calibration.
        The signal will be amplified by a factor
        2**gain before conversion.  A gain switch
        always requires a new calibration.
    updateRate : {AD7705_50HZ, AD7705_60HZ, AD7705_500HZ}
        The input conversion rate for which calibration
        needs to be done. When changing updateRate, the
        ADC always needs a new calibration.
    **/
    AD7705startIO(chipSelectPin);
    SPI.transfer(AD7705_REG_CLOCK | channel);
    SPI.transfer(AD7705_CLK | updateRate);     // Set clock and update rate
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_SELFCAL | gain << 3);  // Start self-calibration
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0xffff);                    // Read data register to ensure DRDY is off.
    AD7705endIO(chipSelectPin);
}

int32_t AD7705getCalibrationRegister(byte chipSelectPin, byte channel,
                                     byte calibrationRegister) {
    /**Return the offset or gain calibration register of an AD7705.
    
    This function can be called after a call to
    AD7705selfCalibrate() to retrieve the calibration
    result, with a short delay(1) in between.
    
    Parameters
    ----------
    chipSelectPin : byte
        The Arduino pin acting as chip-select for the
        ADC from which the register has to be read
    channel : byte
        The ADC channel whose calibration register
        has to be read
    calibrationRegister : {AD7705_REG_OFFSET, AD7705_REG_GAIN}
        The register whose contents are to be read.
        No check is done that in fact calibrationRegister
        is one of the admissible values!                                        // TODO: perhaps we should check this
    */
    AD7705startIO(chipSelectPin);
    SPI.transfer(calibrationRegister | AD7705_READ_REG | channel);
    int32_t result = (int32_t)SPI.transfer16(0xffff) << 8; //read high 16 bits
    result |= SPI.transfer(0xff);              //read low 8 bits
    AD7705endIO(chipSelectPin);
    return result;
}

/** Sets the value in an AD7705 offset or gain calibration register.
 *  'theRegister' must be AD7705_REG_OFFSET or AD7705_REG_GAIN.
 *  This function can be called to restore the result of a previous self-calibration */
int32_t AD7705setCalibrationRegister(byte chipSelectPin, byte channel, byte theRegister, int32_t value) {
  AD7705startIO(chipSelectPin);
  SPI.transfer(theRegister | channel);
  SPI.transfer16(value>>8);                  //write high 16 bits
  SPI.transfer(value&0xff);                  //write low 8 bits
  AD7705endIO(chipSelectPin);
}

/** Waits until the AD7705 with the given chip select pin has finished self-calibration */
void AD7705waitForCalibration(byte chipSelectPin, byte channel) {               // TODO: how does this behave if the power is turned off in the meanwhile?
    AD7705startIO(chipSelectPin);
    do {
        delayMicroseconds(10);
        SPI.transfer(AD7705_REG_SETUP | AD7705_READ_REG | channel);
    } while (SPI.transfer(0x0) & AD7705_SELFCAL);  //read setup register and check for self-calibration
    AD7705endIO(chipSelectPin);
}

/** Waits for data and then reads the given channel of the AD7705 with the given chip select pin.
 *  Returns a signed 16-bit int with the signed value (N.B. we use bipolar mode only) */
int16_t AD7705waitAndReadData(byte chipSelectPin, byte channel) {               // TODO: how does this behave if the power is turned off in the meanwhile?
    while (AD7705readCommRegister(chipSelectPin, channel) & AD7705_DRDY) {//DRDY bit is 0 when ready
        delayMicroseconds(100);
    }
    AD7705startIO(chipSelectPin);
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    int result = SPI.transfer16(0xffff);
    result ^= 0x8000; //flip highest bit: AD7705 output is 0 for -Vref, 0x8000 for 0, 0xffff for +Vref
    AD7705endIO(chipSelectPin);
    return result;
}

/** Reads and returns the communications register of the AD7705 with given chip select pin.
 *  This operation should be done for the channel currently in use. */
byte AD7705readCommRegister(byte chipSelectPin, byte channel) {
    AD7705startIO(chipSelectPin);
    //SPI.transfer16(0xffff); //reset should not be necessary
    //SPI.transfer16(0xffff);
    SPI.transfer(AD7705_READ_REG | AD7705_REG_COMM | channel);
    byte result = SPI.transfer(0xff);
    AD7705endIO(chipSelectPin);
    return result;
}

/* ---------- DAC AD5683 FUNCTIONS ---------- */

/** Resets the AD5683 DAC */
void AD5683reset(byte chipSelectPin) {
    SPI.beginTransaction(AD5683_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
    SPI.transfer(AD5683_RESET_MSB);
    SPI.transfer16(0);
    digitalWrite(chipSelectPin, HIGH);
}

/** Sets the AD5683 DAC output voltage to an unsigned 16-bit value */
void AD5683setVoltage(byte chipSelectPin, uint16_t dacValue) {
    uint16_t lowBytes = dacValue << 4;
    byte hiByte = (dacValue >> 12) | AD5683_SET_DAC;
    SPI.beginTransaction(AD5683_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
    SPI.transfer(hiByte);
    SPI.transfer16(lowBytes);
    digitalWrite(chipSelectPin, HIGH);
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
