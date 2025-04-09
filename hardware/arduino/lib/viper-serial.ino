/* ViPErLEED - Serial communication functions for hardware controllers

@author: Bernhard Mayr
@author: Michael Schmid (@schmid-iap)
@author: Michele Riva (@michele-riva)
@author: Florian Dörr (@FlorianDoerr)
@author: Tun Sinner (@SinTu404)
Date: 26.03.2025
*/

// Libraries
#include <Arduino.h>
#include <EEPROM.h>
#include <SPI.h>
#include <stdarg.h>

/**
This module contains functions for serial communication between PC and
Arduino.

Note that the "raise" function only sets currentState to STATE_ERROR
if STATE_ERROR has been defined. The "Goes to state" documentation
therefore only applies to state machines that use this module and that
contain the STATE_ERROR. currentState and STATE_ERROR are defined in the
states-def.h module.
**/

/** ----------------------------- COMMUNICATION ---------------------------- **/

void raise(byte error_code){
    /**Bring the system to a STATE_ERROR with a given error_code.

    Parameters
    ----------
    error_code
        Byte that identifies the error, see header.

    Writes
    ------
    currentState

    Goes to state
    -------------
    STATE_ERROR : error_code
    **/
#ifdef STATE_ERROR
    errorTraceback[0] = currentState;
    currentState = STATE_ERROR;
#endif // Skip setting state to STATE_ERROR if no state machine is present.
    errorTraceback[1] = error_code;
}


void encodeAndSend(byte *byteArray, uint16_t numBytesBeforeEncoding){
/*
 * Prepares message before sending it to the PC. Changes every
 * byte which happens to have the same value as a MSG_START, and MSG_END or
 * a MSG_SPECIAL_BYTE to two bytes with a leading "MSG_SPECIAL_BYTE" and
 * a following "byte - MSG_SPECIAL_BYTE."
 *
 * Parameters:
 * -----------
 * byteArray : byte*
 *     Pointer to message to be sent
 */
    if (numBytesBeforeEncoding >= MSG_SPECIAL_BYTE){
        raise(ERROR_MSG_SENT_TOO_LONG);
        return;
        }
    byte encodedMessage[2*numBytesBeforeEncoding]; // Worst-case: each byte encoded as two
    byte numBytesAfterEncoding = 0;
    for(int i=0; i < numBytesBeforeEncoding; i++){
        if (byteArray[i] >= MSG_SPECIAL_BYTE){
            encodedMessage[numBytesAfterEncoding] = MSG_SPECIAL_BYTE;
            numBytesAfterEncoding++;
            encodedMessage[numBytesAfterEncoding] = byteArray[i] - MSG_SPECIAL_BYTE;
        }
        else {
            encodedMessage[numBytesAfterEncoding] = byteArray[i];
        }
      numBytesAfterEncoding++;
    }
/*  Send byte array "encodedMessage" (i.e., the actual message) to PC as:
       * MSG_START
       * numbers of bytes in actual message (before encoding, excl. itself and markers)
       * actual message
       * MSG_END
 */
    Serial.write(MSG_START);
    Serial.write(numBytesBeforeEncoding);
    Serial.write(encodedMessage, numBytesAfterEncoding);
    Serial.write(MSG_END);
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
 * Parameters
 * ----------
 * singleByte : byte
 *     The one byte to be sent
 */
  byte byteArray[] = {singleByte};
  encodeAndSend(byteArray, 1);
}


bool checkIfTimedOut(){
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
        raise(ERROR_TIMEOUT);
        return true;
        }
    return false;
}


void debugMsg(const char *message, ...){  // can be a format string
    /** Send a debug message to the PC.

    Caution: when sending a debugMsg to the PC, it is a good idea
    to do so BEFORE sending out a PC_OK, as the python side may be
    triggered into doing some other stuff (including sending other
    serial messages) when a PC_OK is received. Sending a debugMsg
    can disrupt the progress or hold the serial line busy, preventing
    new messages to be actually sent without this being noticed.

    Parameters
    ----------
    message:
        string to be sent. Can contain format characters
    ...:
        variable number of arguments, interpreted as the values
        to be formatted into message.

    First a PC_DEBUG is sent to the PC, then the actual message is
    sent. The message is formatted like "<message % ...>\0", and
    should be at most MSG_SPECIAL_BYTE - 1 characters long, including
    the terminating \0. It is encoded like all others.
    **/
    va_list args;
    va_start(args, message);

    uint16_t n_chars;
    // Note that the buffer size is larger than the maximal message
    // length of MSG_SPECIAL_BYTE - 1
    char _buffer[255];

    n_chars = vsnprintf(_buffer, 255, message, args);
    va_end(args);

    encodeAndSend(PC_DEBUG);
    encodeAndSend(reinterpret_cast<byte*>(_buffer), n_chars);
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
    data_received, msgLength

    Returns
    -------
    True if the message is acceptable

    Goes to state
    -------------
    (unchanged)
        if message is acceptable
    STATE_ERROR : ERROR_MSG_INCONSISTENT
        if the number of decoded bytes does not match the
        length expected from the value that came with the
        message itself
    STATE_ERROR : ERROR_MSG_UNKNOWN
        if the message is an unknown command, or if we
        got some 'data' while we were not expecting any
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        if the message only contains a length byte with
        length 0
    */
    // Decode the message, starting at the second byte,
    // and going up to numBytesRead - 2 (included). This
    // skips MSG_START, the length, and MSG_END
    byte numDecodedBytes = 0;
    for (byte nthByte = 2; nthByte <= numBytesRead - 2; nthByte++) {
        byte decodedByte = serialInputBuffer[nthByte];
        if (decodedByte == MSG_SPECIAL_BYTE) {
            // The actual character is MSG_SPECIAL_BYTE + the next byte
            nthByte++;
            decodedByte += serialInputBuffer[nthByte];
        }
        data_received[numDecodedBytes] = decodedByte;
        numDecodedBytes++;
    }

    // Check if data was received
    msgLength = serialInputBuffer[1];
    if (msgLength == 0){
      raise(ERROR_MSG_DATA_INVALID);
      return false;
    }

    // Check that the number of bytes decoded fits
    if (msgLength != numDecodedBytes){
        raise(ERROR_MSG_INCONSISTENT);
        return false;
        }

    if (numDecodedBytes > 1){
        // Message is some data
        if (not waitingForDataFromPC){
            // But we're not expecting any
            raise(ERROR_MSG_UNKNOWN);
            return false;
        }
        return true;
        // Defer checking to state handlers
    }

	return isAllowedCommand();
}


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
    numBytesRead, serialInputBuffer, readingFromSerial, msgLength,
	newMessage

    Goes to state
    -------------
    STATE_ERROR : ERROR_SERIAL_OVERFLOW
        In case the Arduino serial buffer reaches its limit
    STATE_ERROR : ERROR_MSG_INCONSISTENT
        In case the number of bytes read does not fit with the
        number expected from the first byte after MSG_START

    */
    // Do something only if there is data on the serial line
    if(not Serial.available())  return;

    if(Serial.available() >= SERIAL_BUFFER_SIZE){  // Should never be '>', but better safe than sorry
        // The serial buffer is full and it potentially
        // already discarded some of the bytes that came.
        // The buffer will be flushed in the error handler
        raise(ERROR_SERIAL_OVERFLOW);
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
                raise(ERROR_MSG_TOO_LONG);
                return;
            }
            serialInputBuffer[numBytesRead] = byteRead;
            numBytesRead++;
        }

        // A full message has been read
        if (byteRead == MSG_END) {
            readingFromSerial = false;
            /*for (int i; i < numBytesRead; i++)  // echo message, for debug
                Serial.write(serialInputBuffer[i]);*/
            newMessage = decodeAndCheckMessage();
            return;
        }

        // Delay a very little bit to make sure new characters
        // come in, if any. This should be enough to read in a
        // whole message. If it isn't the case, we will finish
        // during the next state-loop iteration anyway.
        delayMicroseconds(20);
    }
}


/** ---------------------- Serial number methods --------------------- **/

void writeSerialNR(byte *serial_nr) {
	/**Writes the assigned serial number to the EEPROM.**/
    // Check if address only contains allowed values.
    int address = 0;
    byte tmp_char;
    while(address <= 3){
      tmp_char = serial_nr[address];
      if (not ((tmp_char > 47 and tmp_char < 58)
                || (tmp_char > 64 and tmp_char < 91))){
        raise(ERROR_MSG_DATA_INVALID);
        return;
      }
      address += 1;
    }
	// Serial number is stored on EEPROM bytes with addresses 0 to 3.
    address = 0;
    while(address <= 3){
      EEPROM.update(address, serial_nr[address]);
      address += 1;
    }
}


void getSerialNR(byte *serial_nr) {
	/**Get the serial number assigned to the controller.**/
	int address = 0;
    while(address <= 3){
      serial_nr[address] = EEPROM.read(address);
      address += 1;
    }
}
