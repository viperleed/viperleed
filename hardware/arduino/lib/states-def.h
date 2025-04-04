/* ViPErLEED - State machine definitions for hardware controllers

@author: Bernhard Mayr
@author: Michael Schmid (@schmid-iap)
@author: Michele Riva (@michele-riva)
@author: Florian Dörr (@FlorianDoerr)
@author: Tun Sinner (@SinTu404)
Date: 31.03.2025
*/

#ifndef _VIPERLEED_STATES
#define _VIPERLEED_STATES

/** ------------------------- Finite state machine ------------------------- **/
/* In order to create a state machine, users must add their own states
with the corresponding byte value for identification. This can be
achieved in the format:
#define STATE_*   byte_value
States should not be added in this file, but in the .h of the created
state machine.
*/
#define STATE_IDLE                 0  // Wait for requests from PC
#define STATE_ERROR                9  // An error occurred
uint16_t currentState = STATE_IDLE;   // Keeps track of the current state

/** ------------------------- Communication with PC ------------------------ **/
/* For communication purposes, users must define commands and error
codes. Commands are supposed to be single byte messages that instruct
the state machine to enter a certain state or perform a specific task.
Error codes are single byte definitions, that are meant to carry
information about which error occurred. These can be defined in the
following manner:
#define ERROR_*   byte_value (error)
#define PC_*   byte_value (coomand)
Commands and error codes should not be added in this file, but in the
.h of the created state machine.
*/
// Error codes
#define ERROR_RUNTIME           255 /* Some function has been called from an
inappropriate state. This is to flag possible bugs for future development.*/


#endif // _VIPERLEED_STATES
