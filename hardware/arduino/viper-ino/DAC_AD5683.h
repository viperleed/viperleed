/*
ViPErLEED - Driver for AD56836 Digital-to-Analog Converter
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 26.04.2021
---------------------
*/


#ifndef _VIPERLEED_AD5683
#define _VIPERLEED_AD5683

#ifndef _VIPERLEED
    #error Must be included after "viper-ino.h"
#endif

#include <SPI.h>

// Definitions for AD5683 digital-to-analog converter
#define AD5683_SPIMODE   1    // Data accepted at falling edge of SCLK
#define AD5683_SET_DAC   0x30 // Bits of highest (first) byte for writing & setting DAC
                              // Lower 4 bits must be highest 4 bits of data
#define AD5683_RESET_MSB 0x40 // Bits of highest (first) byte for reset (internal reference)


// SPI communication settings
SPISettings AD5683_SPI_SETTING(2500000, MSBFIRST, AD5683_SPIMODE);


#endif