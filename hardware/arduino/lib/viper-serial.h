/* ViPErLEED - Serial communication definitions for hardware controllers

@author: Bernhard Mayr
@author: Michael Schmid (@schmid-iap)
@author: Michele Riva (@michele-riva)
@author: Florian Dörr (@FlorianDoerr)
@author: Tun Sinner (@SinTu404)
Date: 26.03.2025
*/
#ifndef _VIPERLEED_SERIAL
#define _VIPERLEED_SERIAL

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
#define PC_DEBUG            252  // Header for debug messages sent to the PC. No debug can come from the PC
#define PC_ERROR            253  // An error occurred
#define PC_OK                75  // Acknowledge request from PC (ASCCI 'K')

// Error codes
#define ERROR_NO_ERROR            0   // No error
#define ERROR_SERIAL_OVERFLOW     1   // Hardware overflow of Arduino serial
#define ERROR_MSG_TOO_LONG        2   // Too many characters in message from PC
#define ERROR_MSG_INCONSISTENT    3   // Message received from PC is inconsistent. Probably corrupt.
#define ERROR_MSG_UNKNOWN         4   // Unknown request from PC
#define ERROR_MSG_DATA_INVALID    5   // Request from PC contains invalid information
#define ERROR_TIMEOUT             7   // Timed out while waiting for something
#define ERROR_MSG_SENT_TOO_LONG 254   // Any message must be shorter than MSG_SPECIAL_BYTE
byte errorTraceback[2];               // Keeps track of: (0) the state that produced the error, (1) which error occurred (one of ERROR_*)

// Variables used while communicating with the PC
boolean newMessage = false;              // True when a complete, acceptable message has been read
boolean readingFromSerial = false;       // True while Arduino is receiving a message
boolean waitingForDataFromPC = false;       // Keeps track of whether we are waiting for the PC to send something
byte data_received[MSG_MAX_LENGTH];      // Real message. NB: it will never reach MSG_MAX_LENGTH because of markers, length, and encoding
byte msgLength = 0;                      // No. bytes to be expected in message (2nd byte of received message). Used in states that expect data to check the message
byte numBytesRead = 0;                   // Counter for no. of bytes received. numBytesRead may be larger than the number of true data bytes due to encoding.
byte serialInputBuffer[MSG_MAX_LENGTH];  // Contains all the raw (i.e., still encoded) bytes read from the serial line


/** -------------------- Globals for firmware functions -------------------- **/
// Timers (defined in milliseconds)
#define TIMEOUT 4000                  // Max 4 seconds to do stuff
unsigned long initialTime;            // System time when switching to a new state

bool isAllowedCommand(); /** Check if the received 1-byte command is allowed.
	This function needs to reimplemented in state machines to check
	whether any of the received commands is allowed. If an allowed
	command has been received, then the function needs to return true.
	Otherwise, the default return must be false. Similar to:
    switch(data_received[0]){
		{Add allowed commands here in reimplementations.}
        default:
            raise(ERROR_MSG_UNKNOWN);
            return false;
	}
	return true; **/

#include "viper-serial.ino"

#endif //_VIPERLEED_SERIAL
