/*
ViPErLEED - Firmware for Arduino hardware controller
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 31.03.2025
---------------------
*/

#ifndef _VIPERLEED_STATES
#define _VIPERLEED_STATES

/** ------------------------- Finite state machine ------------------------- **/
#define STATE_IDLE                 0  // Wait for requests from PC
#define STATE_ERROR                9  // An error occurred
uint16_t currentState = STATE_IDLE;   // Keeps track of the current state


#endif // _VIPERLEED_STATES
