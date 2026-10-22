/*
ViPErLEED - Driver for AD7705 Analog-to-Digital Converter
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 26.04.2021
---------------------
*/


#ifndef _VIPERLEED_AD7705
#define _VIPERLEED_AD7705

#ifndef _VIPERLEED
    #error Must be included after "viper-ino.h"
#endif



#include <SPI.h>

#define AD7705_DELAY_MICRO 100

// Definitions for AD7705 analog-to-digital converter
#define AD7705_SPIMODE     3   // Data accepted at rising edge of SCLK
#define AD7705_MAX_GAIN    7   // Gains are 0 (x1) to 7 (x128)

//AD7705 communication register
#define AD7705_READ_REG   0x08  // Set bit for reading a register
#define AD7705_CH0        0x00  // Set bit for addressing channel 0
#define AD7705_CH1        0x01  // Set bit for addressing channel 1
#define AD7705_REG_COMM   0x00  // Bits to set for accessing communication register
#define AD7705_REG_SETUP  0x10  // Bits to set for accessing setup register
#define AD7705_REG_CLOCK  0x20  // Bits to set for accessing clock register
#define AD7705_REG_DATA   0x30  // Bits to set for accessing data register (16 bit)
#define AD7705_REG_OFFSET 0x60  // Bits to set for accessing calibration offset register (24 bit)
#define AD7705_REG_GAIN   0x70  // Bits to set for accessing calibration gain register (24 bit)
#define AD7705_DRDY       0x80  // Bit mask for Data Ready bit in communication register

//AD7705 setup register
#define AD7705_TRIGGER    0x01  // Bit to set for triggering ADC (called FSYNC in data sheet)
#define AD7705_BUFMODE    0x02  // Bit to set for buffered mode (not used by us)
#define AD7705_UNIPOLAR   0x04  // Bit to set for unipolar mode (not used by us)
#define AD7705_SELFCAL    0x40  // Bit to set for self-calibration

//AD7705 clock register
#define AD7705_CLK        0x0C  // Bits to set for 4.9152 MHz crystal: CLK=1, CLKDIV=1
#define AD7705_50HZ       0x04  // Bits to set for 50 Hz update rate and suppression
#define AD7705_60HZ       0x05  // Bits to set for 60 Hz update rate and suppression
#define AD7705_500HZ      0x07  // Bits to set for 500 Hz update rate

// Correction for the finite source impedance and ADC input impedance
#define R_IN_0              4.1e6       // input resistance of AD7705 in gain0 = x1, measured typical value
#define REF_OVER_RANGE  (2500.0/32768)  // millivolts (uA@1kOhm) per bit at gain0, bipolar, input 0ohm
#define R_SOURCE         1300.0         // Source resistance of our circuit at ADC inputs, 1.3 kOhm.

const float millivoltsPerBit[] = {
    REF_OVER_RANGE/R_IN_0*(R_IN_0 + R_SOURCE),         // gain0 = x1
    REF_OVER_RANGE/R_IN_0*2*(R_IN_0/2 + R_SOURCE)/2,   // gain1 = x2 has half R_IN_0
    REF_OVER_RANGE/R_IN_0*4*(R_IN_0/4 + R_SOURCE)/4,   // gain2 = x4 has 1/4 R_IN_0
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/8,   // gain3 = x8 and up: 1/8 R_IN_0
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/16,  // gain4 = x16
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/32,  // gain5 = x32
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/64,  // gain6 = x64
    REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/128, // gain7 = x128
};


// SPI communication settings
SPISettings AD7705_SPI_SETTING(2500000, MSBFIRST, AD7705_SPIMODE);







#endif
