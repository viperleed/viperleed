/*
ViPErLEED - Driver for TLE7209-3R H-Bridge
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 21.04.2023
---------------------
*/


#ifndef _VIPERLEED_TLE7209
#define _VIPERLEED_TLE7209


#include <SPI.h>


#define TLE7209_DELAY_MICRO 100

/*
Communication parameters for the TLE7209-3R:
  - The TLE7209 always operates in slave mode (read-only)
  - Baud rate: 2 MBaud/s max.
  - MSbit first, clock polarity (CPOL=0) and phase (CPHA=1)
*/
#define TLE7209_SPI_BAUD  1E5        // Can be increased later, up to 2E6
#define TLE7209_SPIMODE   SPI_MODE1  // SPI_MODE1 means CPOL=0, CPHA=1
#define TLE7209_SPI_SETTINGS SPISettings(TLE7209_SPI_BAUD, MSBFIRST, TLE7209_SPIMODE)


// TLE7209 SPI instruction-byte encoding
#define TLE7209_READ_IDENTIFIER     0x00
#define TLE7209_READ_VERSION        0x03
#define TLE7209_READ_DIAG_REGISTER  0x09
#define TLE7209_SPI_TRANSMISSION_OK 0b00101010
#define TLE7209_DEFAULT_DEVICE_ID   0b10100010


// TLE7209 diagnostic register definitions
#define TLE7209_EN_DIS    0x80	// EN/DIS = 0 if EN = low or DIS = high
#define TLE7209_OVER_TEMP 0x40	// OT = 0 in case of over-temperature
#define TLE7209_CURR_RED  0x20	// CurrRed = 0 in case of temperature-dependent current limitation
#define TLE7209_CURR_LIM  0x10	// CurrLim = 0 in case of current limitation
#define TLE7209_DIA_21 0x08	// Diagnostic bit 2 of output OUT2
#define TLE7209_DIA_20 0x04	// Diagnostic bit 1 of output OUT2
#define TLE7209_DIA_11 0x02	// Diagnostic bit 2 of output OUT1
#define TLE7209_DIA_10 0x01	// Diagnostic bit 1 of output OUT1
#define TLE7209_ALL_ERROR_BITS 0b01111111 // All bits above, except the MSB


// TLE7209 error codes
enum TLE7209_Error {
    TLE7209_NoError,
    TLE7209_TransmissionError,
    TLE7209_InvalidDeviceId,
    TLE7209_DiagnosticsError,
};


// TLE7209 verification byte definitions
#define TLE7209_TRANS_F 0x01	// Bit is set if previous transfer was recognized as valid
                                // Bit is cleared on error during previous transfer


#endif
