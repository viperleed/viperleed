/*
ViPErLEED - Firmware for Arduino hardware controller
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 26.04.2021
---------------------
*/

#ifndef _VIPERLEED
#define _VIPERLEED

#include "ADC_AD7705.h"  // Settings of the currently used ADC
#include "DAC_AD5683.h"  // Settings of the currently used DAC


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
#define MSG_SPECIAL_BYTE 252       // Prevents clashing of data with (start, end, error)
#define MSG_MAX_LENGTH 16          // Max no. of bytes in an incoming serial message     // TODO: NEED TO INCREASE THIS TO HOST THE 3 FLOAT ADC VALUES
#ifndef SERIAL_BUFFER_SIZE
    #define SERIAL_BUFFER_SIZE 64  // Arduino limit for serial buffer. If buffer is full, new bytes are DISCARDED
#endif

// Acceptable messages for communication with the PC
#define PC_AUTOGAIN       65    // PC requested auto-gain for ADCs (ASCII 'A')
#define PC_CALIBRATION    67    // PC requested self-calibration of all ADCs at all gains (ASCII 'C')
#define PC_ERROR         253    // An error occurred
#define PC_CONFIGURATION  63    // PC requested hardware configuration (ASCII '?')
#define PC_SET_UP_ADCS    83    // PC requested to prepare the ADCs for a measurement (ASCII 'S')   // TODO: requires calibration data up to date (keep a flag, raise errors if never calibrated since bootup. Flag should be set to false at reset(). The python side will have to keep track of when the last calibration was done, and warn if it is too old.)
#define PC_OK             75    // Acknowledge request from PC (ASCCI 'K')
#define PC_RESET          82    // PC requested a global reset (ASCII 'R')
#define PC_SET_VOLTAGE    86    // PC requested to set a certain energy (ASCII 'V')

// Error codes
#define ERROR_NO_ERROR            0   // No error
#define ERROR_SERIAL_OVERFLOW     1   // Hardware overflow of Arduino serial
#define ERROR_MSG_TOO_LONG        2   // Too many characters in message from PC
#define ERROR_MSG_INCONSITENT     3   // Message received from PC is inconsistent. Probably corrupt.
#define ERROR_MSG_UNKNOWN         4   // Unknown request from PC
#define ERROR_MSG_DATA_INVALID    5   // Request from PC contains invalid information
#define ERROR_NEVER_CALIBRATED    6   // The ADCs have never been calibrated before since bootup
#define ERROR_TIMEOUT             7   // Timed out while waiting for something
#define ERROR_ADC_SATURATED       8   // One of the ADC values reached saturation, and gain can't be decreased further
#define ERROR_TOO_HOT             9   // The temperature read by the LM35 is too high
#define ERROR_RUNTIME           255   // Some function has been called from an inappropriate state. This is to flag possible bugs for future development.
byte errorTraceback[2];               // Keeps track of: (0) the state that produced the error, (1) which error occurred (one of ERROR_*)

// Variables used while communicating with the PC
byte numBytesRead = 0;                   // Counter for no. of bytes received. numBytesRead may be larger than the number of true data bytes due to encoding.
byte msgLength = 0;                      // No. bytes to be expected in message (2nd byte of received message). Used in states that expect data to check the message
byte serialInputBuffer[MSG_MAX_LENGTH];  // Contains all the raw (i.e., still encoded) bytes read from the serial line
boolean readingFromSerial = false;       // True while Arduino is receiving a message
boolean newMessage = false;              // True when a complete, acceptable message has been read

/** TODO:
    Judging from the code, and if I didn't misinterpret anything, I think
    we can have only one buffer for the input/output to the PC, so we can
    replace "data_received" and "data_send" with a single "message" array
    Alternatively, we can have them be a decodedCommand[] and encodedReply[]
*/
byte data_received[MSG_MAX_LENGTH];      // Real message
byte data_send[MSG_MAX_LENGTH];




/** ------------------------- Finite state machine ------------------------- **/

#define STATE_IDLE                 0  // Wait for requests from PC
#define STATE_SET_UP_ADCS          1  // Pick correct ADC channels and no. of measurement points
#define STATE_SET_VOLTAGE          2  // Set a voltage with the DAC and triggers the ADCs
#define STATE_MEASURE_ADCS         4  // ADC measurements in progress
#define STATE_ADC_VALUES_READY     5  // ADC measurements done
#define STATE_AUTOGAIN_ADCS        6  // Find optimal gain for both ADCs
#define STATE_GET_CONFIGURATION    7  // Find current hardware configuration and return it with the firmware version
#define STATE_CALIBRATE_ADCS       8  // Figure out correct offset and calibration factors for ADCs at all gains.
#define STATE_ERROR                9  // An error occurred
uint16_t currentState = STATE_IDLE;   // Keeps track of the current state
bool waitingForDataFromPC = false;    // Keeps track of whether we are in a state that is waiting for the PC to send something




/** ---------------------- Hardware-specific settings ---------------------- **/
// Measurement devices
#define N_MAX_MEAS         3  // Number of measurement devices that we can have: ADC#0, ADC#1, LM35 (if present)
#define N_MAX_ADCS_ON_PCB  2  // Maximum number of ADCs on the PCB

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



/** -------------------- Globals for firmware functions -------------------- **/
// Timers (defined in milliseconds)
#define TIMEOUT 5000                  // Max 5 seconds to do stuff
unsigned long initialTime;            // System time when switching to a new state
uint16_t      dacSettlingTime = 100;  // The time interval for the DAC output to be stable
                                      //   This is just a default value. The actual one comes
                                      //   from the PC with a PC_SET_VOLTAGE command, and is
                                      //   read in setVoltage()

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
byte     calibrationGain = 0;                         // Gain for which a calibration is currently being performed (in parallel for both ADCs)
bool     calibratedChannels[N_MAX_ADCS_ON_PCB][2];    // Array of flags to keep track of which channels (last index) of the ADCs have been calibrated
int32_t  selfCalDataVsGain[AD7705_MAX_GAIN + 1][2][2][2]; // For each gain, two ADCs, two channels each, and last index is offset(0)&gain(1)

// ADCs: variables for measuring and storing the line frequency ripple
int16_t  maximumPeak[N_MAX_ADCS_ON_PCB];       // Maximum of measurement, one for each ADC
int16_t  minimumPeak[N_MAX_ADCS_ON_PCB];       // Minimum of measurement, one for each ADC
int16_t  adc0RipplePP = 0;     // Ripple (peak-peak) measured at gain=0 for ADC#0 during auto-gain  // TODO: this is unused but should be
int16_t  adc1RipplePP = 0;     // Ripple (peak-peak) measured at gain=0 for ADC#1 during auto-gain  // TODO: this is unused but should be

/* // ADC container                                                             // TODO: use it to simplify code below
struct analogToDigitalConverter {
    byte     chipSelect;    // to be initialized!
    uint16_t present;       // to be initialized!
    byte     channel = AD7705_CH0;
    byte     gain = 0;
    bool     shouldDecreaseGain = false;
    int16_t  ripplePP = 0;                                   // Ripple (peak-peak) measured at gain=0 during auto-gain    // TODO: should we keep the largest ripple value ever measured or the latest?
    int32_t  calibrationForGain[AD7705_MAX_GAIN + 1][2][2];  // For each gain, two channels each, and last index is offset(0) & gain(1)
	bool     calibratedChannels[2] = {false, false};
} externalADCs[N_MAX_ADCS_ON_PCB]; */

// Measurements
uint16_t numMeasurementsToDo = 1;         // No. of ADC measurements to do before measurement is considered over
uint16_t numMeasurementsToDoBackup = 1;   // Copy of the previous one, used to restore the previous value after auto-gain is done
uint16_t numMeasurementsDone;             // Counter for the number of ADC measurements done so far
int32_t  summedMeasurements[N_MAX_MEAS];  // Measurements of ADCs and LM35 are summed up here

floatOrBytes fDataOutput[N_MAX_MEAS];     // Measurements in physical units  // TODO: rename







#endif
