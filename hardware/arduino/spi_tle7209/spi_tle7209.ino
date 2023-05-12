/*
Test firmware for Arduino SPI driver module (TLE7209.ino/TLE7209.h)
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 12.05.2023
---------------------
*/

// Libraries
#include <SPI.h>
#include "spi_tle7209.h"      // Arduino-related settings, also includes motor driver header

#define DEBUG true            // Debug mode, writes to serial line, for use in serial monitor
#define TLE_CHIPSELECT 11     // PB7 on the Arduino Micro board


void setup() {
    #if DEBUG
        Serial.setTimeout(100);
        Serial.begin(9600);   // opens serial port, sets data rate to 9600 bps
    #endif

    pinMode(TLE_CHIPSELECT, OUTPUT);
    SPI.begin();              // Initializes the SPI bus (SCK and MOSI as OUTPUT)
    pinMode(MISO, INPUT);     // MISO = pin PB3
}


void loop() {
    uint8_t byteRead;
    TLE7209_Error errcode = TLE7209_NoError;
    
    #if DEBUG
        Serial.println("Heartbeat\n");
    #endif

    errcode = TLE7209readIDandVersion(TLE_CHIPSELECT, &byteRead);
    delayMicroseconds(50);

    byteRead = TLE7209readDiagnosticRegister(TLE_CHIPSELECT, &byteRead);
    delayMicroseconds(50);
}

