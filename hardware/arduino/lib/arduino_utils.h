/** Useful functions for Arduino, not included in STDlib

@author: Michele Riva (@michele-riva)
@author: Michael Schmid (@schmid-iap)
@author: Florian Dörr (@FlorianDoerr)
@author: Christoph Pfungen (@cpfungen)
@author: Stefan Mitterhöfer (@Stefan-Mitterhoefer)
@author: Tun Sinner (@SinTu404)
**/

#ifndef _VIPERLEED_ARDUINO_UTILS
#define _VIPERLEED_ARDUINO_UTILS


#include <Arduino.h>  // For pin definitions and macros
#include <SPI.h>


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

uint16_t bigger(uint16_t a, uint16_t b) {
    return (a > b) ? a : b;
}

int32_t bigger(int32_t a, int32_t b) {
    return (a > b) ? a : b;
}

uint16_t biggest(uint16_t a, uint16_t b, uint16_t c) {
    return bigger(a, bigger(b, c));
}

int32_t biggest(int32_t a, int32_t b, int32_t c) {
    return bigger(a, bigger(b, c));
}

/** Gets the median of three numbers */
uint16_t getMedian(uint16_t a0, uint16_t a1, uint16_t a2) {
    uint16_t maximum = biggest(a0, a1, a2);
    if (maximum == a0) return bigger(a1, a2);
    if (maximum == a1) return bigger(a0, a2);
    else return bigger(a0, a1);
}

/** Gets the median of three numbers */
int32_t getMedian(int32_t a0, int32_t a1, int32_t a2) {
    int32_t maximum = biggest(a0, a1, a2);
    if (maximum == a0) return bigger(a1, a2);
    if (maximum == a1) return bigger(a0, a2);
    else return bigger(a0, a1);
}

uint16_t analogReadMedian(byte pin) {
    /**
    Read the Arduino ADC for a given pin
    (A0...) three times and return the median
    */
    uint16_t a0 = analogRead(pin);
    uint16_t a1 = analogRead(pin);
    uint16_t a2 = analogRead(pin);
    return getMedian(a0, a1, a2);
}

#define sgn(x) (x < 0 ? -1 : 1)

#define log2(x) (log(x)/log(2))

#endif
