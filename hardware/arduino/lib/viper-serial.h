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

/** ------------------------ Communication with PC ----------------------- **/
// Constants for communication with the PC

// MSG_START: Beginning of a serial message
#define MSG_START             254
// MSG_END: End of a serial message
#define MSG_END               255
// MSG_SPECIAL_BYTE: Prevents clashing of data with
// MSG_START, MSG_END, PC_DEBUG, and PC_ERROR.
#ifndef MSG_SPECIAL_BYTE
    #define MSG_SPECIAL_BYTE  251
#endif
// MSG_MAX_LENGTH: Max no. of bytes in an incoming
// serial message. One less than SERIAL_BUFFER_SIZE.
#define MSG_MAX_LENGTH         63
// SERIAL_BUFFER_SIZE: Arduino limit for serial buffer.
// If buffer is full, new bytes are DISCARDED.
#ifndef SERIAL_BUFFER_SIZE
    #define SERIAL_BUFFER_SIZE 64
#endif

// Acceptable messages for communication with the PC

// PC_DEBUG: Header for debug messages sent to the PC.
// No debug can come from the PC.
#define PC_DEBUG            252
// PC_ERROR: An error occurred
#define PC_ERROR            253
// PC_OK: Acknowledge request from PC (ASCCI 'K')
#define PC_OK                75

// Serial-communication-related error codes

// ERROR_NO_ERROR: No error.
#define ERROR_NO_ERROR            0
// ERROR_SERIAL_OVERFLOW: Hardware overflow of Arduino serial.
#define ERROR_SERIAL_OVERFLOW     1
// ERROR_MSG_TOO_LONG: Too many characters in message from PC.
#define ERROR_MSG_TOO_LONG        2
// ERROR_MSG_INCONSISTENT: Message received from
// PC is inconsistent. Probably corrupt.
#define ERROR_MSG_INCONSISTENT    3
// ERROR_MSG_UNKNOWN: Unknown request from PC.
#define ERROR_MSG_UNKNOWN         4
// ERROR_MSG_DATA_INVALID: Request from PC contains invalid information.
#define ERROR_MSG_DATA_INVALID    5
// ERROR_TIMEOUT: Timed out while waiting for something.
#define ERROR_TIMEOUT             7
// ERROR_MSG_SENT_TOO_LONG: Any message must be shorter than MSG_SPECIAL_BYTE.
#define ERROR_MSG_SENT_TOO_LONG 254
// errorTraceback: Keeps track of: (0) the state that produced the error,
// (1) which error occurred (one of ERROR_*).
byte errorTraceback[2];

// Variables used while communicating with the PC

// newMessage: True when a complete, acceptable message has been read.
boolean newMessage = false;
// readingFromSerial: True while Arduino is receiving a message.
boolean readingFromSerial = false;
// waitingForDataFromPC: Keeps track of whether we are
// waiting for the PC to send something.
boolean waitingForDataFromPC = false;
// data_received: Real message. NB: it will never reach
// MSG_MAX_LENGTH because of markers, length, and encoding.
byte data_received[MSG_MAX_LENGTH];
// msgLength: No. bytes to be expected in message (2nd byte of received
// message). Used in states that expect data to check the message.
byte msgLength = 0;
// numBytesRead: Counter for no. of bytes received. numBytesRead may be
// larger than the number of true data bytes due to encoding.
byte numBytesRead = 0;
// serialInputBuffer: Contains all the raw (i.e., still encoded)
// bytes read from the serial line.
byte serialInputBuffer[MSG_MAX_LENGTH];


/** ------------------- Globals for firmware functions ------------------- **/
// Timers (defined in milliseconds)

// TIMEOUT: Max 4 seconds to perform tasks.
#define TIMEOUT 4000
// initialTime: System time, used to set a start point for timeouts.
unsigned long initialTime;

// isAllowedCommand: Check if the received 1-byte command is allowed.
// This function needs to be reimplemented in state machines to check
// whether any of the received commands is allowed. If an allowed
// command has been received, then the function needs to return true.
// Otherwise, the default return must be false. Similar to:
/**
    switch(data_received[0]){
        case <allowed_command_1>: break;
        case <allowed_command_2>: break;
        case <allowed_command_3>: break;
        {Add allowed commands here in reimplementations.}
        default:
            raise(ERROR_MSG_UNKNOWN);
            return false;
    }
    return true; **/
bool isAllowedCommand();

#include "viper-serial.ino"

#endif //_VIPERLEED_SERIAL
