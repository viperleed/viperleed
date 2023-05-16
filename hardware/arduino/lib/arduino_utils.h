/**Useful functions for Arduino, not included in STDlib

@author: Michele Riva
@author: Michael Schmid
@author: Florian Doerr
@author: Christoph Pfungen
**/

#ifndef _VIPERLEED_ARDUINO_UTILS
#define _VIPERLEED_ARDUINO_UTILS


#include <Arduino.h>  // For pin definitions and macros


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


#endif
