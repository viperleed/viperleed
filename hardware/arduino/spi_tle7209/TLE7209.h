/*
ViPErLEED - Driver for TLE7209-3R H-Bridge
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 21.04.2023
---------------------
*/


#ifndef _VIPERLEED_TLE7209
#define _VIPERLEED_TLE7209

// #ifndef _VIPERLEED
//     #error Must be included after "viper-ino.h"
// #endif



#include <SPI.h>

#define TLE7209_DELAY_MICRO 100

// Definitions for TLE7209-3R H-Bridge
#define TLE7209_SPIMODE SPI_MODE1   // SPI_MODE1 should translate to CPOL = 0, CPHA = 1


// TLE7209 SPI instruction-bye encoding
#define TLE7209_READ_IDENTIFIER 0x00
#define TLE7209_READ_VERSION 0x03
#define TLE7209_READ_DIAG_REGISTER 0x09


// TLE7209 diagnostic register definitions
#define TLE7209_EN_DIS   0x80	// EN/DIS = 0 if EN = low or DIS = high 
#define TLE7209_OVER_TEMP 0x40	// OT = 0 in case of over-temperature
#define TLE7209_CURR_RED 0x20	// CurrRed = 0 in case of temperature-dependent current limitation
#define TLE7209_CURR_LIM 0x10	// CurrLim = 0 in case of current limitation
#define TLE7209_DIA21 0x08	// Diagnostic bit 2 of output OUT2
#define TLE7209_DIA_20 0x04	// Diagnostic bit 1 of output OUT2
#define TLE7209_DIA_11 0x02	// Diagnostic bit 2 of output OUT1
#define TLE7209_DIA_10 0x01	// Diagnostic bit 1 of output OUT1
 

// TLE7209 verification byte definitions
#define TLE7209_TRANS_F 0x01	// Bit is set if previous transfer was recognized as valid
				// Bit is cleared on error during previous transfer


// SPI communication settings
SPISettings TLE7209_SPI_SETTING(100000, MSBFIRST, TLE7209_SPIMODE);




#endif
