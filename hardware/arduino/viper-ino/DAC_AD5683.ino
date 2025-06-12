/*
ViPErLEED - Driver for AD56836 Digital-to-Analog Converter
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 26.04.2021
---------------------
*/

#include "DAC_AD5683.h"

/* ---------- DAC AD5683 FUNCTIONS ---------- */

/** Resets the AD5683 DAC */
void AD5683reset(byte chipSelectPin) {
    SPI.beginTransaction(AD5683_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
    SPI.transfer(AD5683_RESET_MSB);
    SPI.transfer16(0);
    digitalWrite(chipSelectPin, HIGH);
}

/** Sets the AD5683 DAC output voltage to an unsigned 16-bit value */
void AD5683setVoltage(byte chipSelectPin, uint16_t dacValue) {
    uint16_t lowBytes = dacValue << 4;
    byte hiByte = (dacValue >> 12) | AD5683_SET_DAC;
    SPI.beginTransaction(AD5683_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
    SPI.transfer(hiByte);
    SPI.transfer16(lowBytes);
    digitalWrite(chipSelectPin, HIGH);
}