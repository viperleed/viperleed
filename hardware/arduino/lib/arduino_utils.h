/**Useful functions for Arduino, not included in STDlib

@author: Michele Riva
@author: Michael Schmid
@author: Florian Doerr
@author: Christoph Pfungen
@author: Stefan Mitterhöfer
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

void SPI_initialize() {
    /**
    Initialize the SPI bus.
    Set Serial Clock SCK (Pin PB1) and MOSI (Pin PB2) aus OUTPUT.
    Set MISO (Pin PB3) as INPUT.
    **/
    SPI.begin();
    pinMode(MISO, INPUT);
}


uint16_t bigger16(uint16_t a, uint16_t b) {
  return (a > b) ? a : b;
}

uint16_t biggest16(uint16_t a, uint16_t b, uint16_t c) {
  return bigger16(a, bigger16(b, c));
}

/** Gets the median of three numbers */
uint16_t getMedian16(uint16_t a0, uint16_t a1, uint16_t a2) {
  uint16_t maximum = biggest16(a0, a1, a2);
  if (maximum == a0) return bigger16(a1, a2);
  if (maximum == a1) return bigger16(a0, a2);
  else return bigger16(a0, a1);
}

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

#define sgn(x) (x < 0 ? -1 : 1)
//static inline int8_t sgn(float val) {
//  if (val < 0) return -1;
////  if (val == 0) return 0;
//  return 1;
//}

#define log2(x) (log(x)/log(2))
//float log2(float val)
//{
//   return log(val) / log(2);
//}



#endif
