/*
ViPErLEED - Firmware for Arduino hardware controller
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr, Tun Sinner
Date: 26.03.2025
---------------------
*/

// Note that the required definitions of currentState and STATE_ERROR
// must be present in the .h file of the importing main module.

/** ------------------------- Communication with PC ------------------------ **/

// Constants for communication with the PC
#define MSG_START         254      // Beginning of a serial message
#define MSG_END           255      // End of a serial message
#define MSG_SPECIAL_BYTE  251      // Prevents clashing of data with (start, end, error)
#define MSG_MAX_LENGTH     63      // Max no. of bytes in an incoming serial message; one less than SERIAL_BUFFER_SIZE
#ifndef SERIAL_BUFFER_SIZE
    #define SERIAL_BUFFER_SIZE 64  // Arduino limit for serial buffer. If buffer is full, new bytes are DISCARDED
#endif

// Acceptable messages for communication with the PC
#define PC_AUTOGAIN          65  // PC requested auto-gain for ADCs (ASCII 'A')
#define PC_CALIBRATION       67  // PC requested self-calibration of all ADCs at all gains (ASCII 'C')
#define PC_CHANGE_MEAS_MODE 109  // PC requested a change between continuous and single measurement mode (ASCII 'm')
#define PC_CONFIGURATION     63  // PC requested hardware configuration (ASCII '?')
#define PC_DEBUG            252  // Header for debug messages sent to the PC. No debug can come from the PC
#define PC_ERROR            253  // An error occurred
#define PC_MEASURE_ONLY      77  // PC requested measurement without changing Voltage (ASCII 'M')
#define PC_OK                75  // Acknowledge request from PC (ASCCI 'K')
#define PC_RESET             82  // PC requested a global reset (ASCII 'R')
#define PC_SET_SERIAL_NR    115  // PC requested serial number (ASCII 's')
#define PC_SET_UP_ADCS       83  // PC requested to prepare the ADCs for a measurement (ASCII 'S')   // TODO: The python side will have to keep track of when the last calibration was done, and warn if it is too old.
#define PC_SET_VOLTAGE       86  // PC requested to set a certain energy (ASCII 'V')
#define PC_SET_VOLTAGE_ONLY 118  // PC requested set energy without follow up measurement (ASCII 'v')
#define PC_STOP             120  // PC requested a stop on all activity. Return to idle (ASCII 'x')

// Error codes
#define ERROR_NO_ERROR            0   // No error
#define ERROR_SERIAL_OVERFLOW     1   // Hardware overflow of Arduino serial
#define ERROR_MSG_TOO_LONG        2   // Too many characters in message from PC
#define ERROR_MSG_INCONSISTENT    3   // Message received from PC is inconsistent. Probably corrupt.
#define ERROR_MSG_UNKNOWN         4   // Unknown request from PC
#define ERROR_MSG_DATA_INVALID    5   // Request from PC contains invalid information
#define ERROR_NEVER_CALIBRATED    6   // The ADCs have never been calibrated before since bootup
#define ERROR_TIMEOUT             7   // Timed out while waiting for something
#define ERROR_ADC_SATURATED       8   // One of the ADC values reached saturation, and gain can't be decreased further
#define ERROR_TOO_HOT             9   // The temperature read by the LM35 is too high
#define ERROR_MSG_SENT_TOO_LONG 254   // Any message must be shorter than MSG_SPECIAL_BYTE
#define ERROR_RUNTIME           255   // Some function has been called from an inappropriate state. This is to flag possible bugs for future development.
byte errorTraceback[2];               // Keeps track of: (0) the state that produced the error, (1) which error occurred (one of ERROR_*)

// Variables used while communicating with the PC
boolean newMessage = false;              // True when a complete, acceptable message has been read
boolean readingFromSerial = false;       // True while Arduino is receiving a message
byte data_received[MSG_MAX_LENGTH];      // Real message. NB: it will never reach MSG_MAX_LENGTH because of markers, length, and encoding
byte msgLength = 0;                      // No. bytes to be expected in message (2nd byte of received message). Used in states that expect data to check the message
byte numBytesRead = 0;                   // Counter for no. of bytes received. numBytesRead may be larger than the number of true data bytes due to encoding.
byte serialInputBuffer[MSG_MAX_LENGTH];  // Contains all the raw (i.e., still encoded) bytes read from the serial line

/** ------------------------- Finite state machine ------------------------- **/
bool waitingForDataFromPC = false;    // Keeps track of whether we are in a state that is waiting for the PC to send something

/** -------------------- Globals for firmware functions -------------------- **/
// Timers (defined in milliseconds)
#define TIMEOUT 4000                  // Max 4 seconds to do stuff
unsigned long initialTime;            // System time when switching to a new state
