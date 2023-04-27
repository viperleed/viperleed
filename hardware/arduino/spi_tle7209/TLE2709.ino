/*
ViPErLEED - Driver for TLE7209-3R H-Bridge
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 21.04.2023
---------------------
*/


/* 
Communication parameters for the TLE7209-3R:
  - The TLE7209 always operates in slave mode (read-only)
  - Max. baud rate: 2 MBaud
  - SPI active when 'DMS' > 3.5V
  - Clock polarity and phase: CPOL = 0, CPHA = 1
*/


#include "TLE7209.h"

#define TLE_CHIPSELECT 11       // PB7


/* ---------- TLE7209 FUNCTIONS ---------- */

/** Starts an I/O operation for the TLE7209 with given chip select pin */
void TLE7209startIO(byte chipSelectPin) {
    SPI.beginTransaction(TLE7209_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
}


/** Finishes an I/O operation for the TLE7209 with given chip select pin */
void TLE7209endIO(byte chipSelectPin) {
    digitalWrite(chipSelectPin, HIGH);
}


byte readTLE7209(byte chipSelectPin, byte request){     // request = {TLE7209_READ_IDENTIFIER, TLE7209_READ_VERSION, TLE7209_READ_DIAG_REGISTER}
  TLE7209startIO(chipSelectPin);
  uint16_t bytesRead = SPI.transfer16((uint16_t)(request << 8));
  TLE7209endIO(chipSelectPin);

  // Check verification byte; If LSbit is clear, last transfer was recognized as valid
  if(((uint8_t)(bytesRead >> 8) | 0B00101010) != 0B00101010)
    return -1; 

  // Return only the LSB (data byte)
  return (uint8_t)bytesRead;
}

// Perform two transfers back-to-back (READ_ID, READ_VER)
 byte TLE7209readIDandVersion(byte chipSelectPin){
  byte deviceID = readTLE7209(chipSelectPin, TLE7209_READ_IDENTIFIER); 
  
  if(deviceID == -1)
    return -2;

  // Chip ID is fixed, check if match
  if(deviceID != 0B10100010)
    return -1;
    
  return readTLE7209(chipSelectPin, TLE7209_READ_VERSION);
}


byte TLE7209readDiagnosticRegister(byte chipSelectPin){
  byte diagnosticReg = readTLE7209(chipSelectPin, TLE7209_READ_DIAG_REGISTER);

  if(diagnosticReg == -1)
    return -2;

  switch(diagnosticReg)   // Handle diagnostic bits, return corresponding error code
  {
    case TLE7209_EN_DIS:
      break;
    case TLE7209_OVER_TEMP:
      break;
    case TLE7209_CURR_RED:
      break;
    case TLE7209_CURR_LIM:
      break;
    case TLE7209_DIA21:
      break;
    case TLE7209_DIA_20:
      break;
    case TLE7209_DIA_11:
      break;
    case TLE7209_DIA_10:
      break;

    default:
    {
      // Handle combinations according to fault priority (TLE7209-3R datasheet, p. 16)
    }
  }

  return diagnosticReg;
}


void setup()
{ 
  pinMode(TLE_CHIPSELECT, OUTPUT);  

  SPI.begin();                    // SPI.begin initializes the SPI bus (defines SCK and MOSI as output pins)
                                  // Alternatively, use pinMode() on SCK and MOSI below
  pinMode(SCK, OUTPUT);           // SCK = pin PB1
  pinMode(MOSI, OUTPUT);          // MOSI = pin PB2
  pinMode(MISO, INPUT);           // MISO = pin PB3
} 


void loop()
{ 
  uint8_t byteRead;

  byteRead = TLE7209readIDandVersion(TLE_CHIPSELECT);  
  delayMicroseconds(50);
  
  byteRead = TLE7209readDiagnosticRegister(TLE_CHIPSELECT);
  delayMicroseconds(50);
}
