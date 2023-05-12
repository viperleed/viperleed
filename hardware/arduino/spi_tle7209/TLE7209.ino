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
  - Baud rate: 2 MBaud/s max.
  - MSbit first, clock polarity and phase: CPOL = 0, CPHA = 1
*/


#include "TLE7209.h"

#define DEBUG 1

#define TLE_CHIPSELECT 11       // PB7                                          // TODO: will need to go in the Arduino firmware, not in the library


/* ---------- TLE7209 FUNCTIONS ---------- */

/** Starts an I/O operation for the TLE7209 with given chip-select pin */
void TLE7209startIO(byte chipSelectPin) {
    SPI.beginTransaction(TLE7209_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
}


/** Finishes an I/O operation for the TLE7209 with given chip-select pin */
void TLE7209endIO(byte chipSelectPin) {
    digitalWrite(chipSelectPin, HIGH);
    SPI.endTransaction();
}

                                                                                // TODO: docstring
TLE7209_Error readTLE7209(byte chipSelectPin, byte request, byte *data){        // request = {TLE7209_READ_IDENTIFIER, TLE7209_READ_VERSION, TLE7209_READ_DIAG_REGISTER}
    TLE7209startIO(chipSelectPin);
    uint16_t bytesRead = SPI.transfer16((uint16_t)(request << 8));
    TLE7209endIO(chipSelectPin);

    // Check verification byte
    uint8_t transmit_ok = bytesRead >> 8;
    transmit_ok &= 0b00111111; // The highest two bits are not relevant
    if(transmit_ok != TLE7209_SPI_TRANSMISSION_OK){
        #if DEBUG
            Serial.println("Verification byte: TRANS_F is set and/or wrong bit toggle sequence detected\n");
        #endif
        return TLE7209_TransmissionError;
    }

    // Return only the LSB (data byte)
    *data = (uint8_t)bytesRead;
    return TLE7209_NoError;
}


// Perform two transfers back-to-back (READ_ID, READ_VER)                       // TODO: docstring
TLE7209_Error TLE7209readIDandVersion(byte chipSelectPin, byte *version){
    byte deviceID;
    TLE7209_Error errcode = TLE7209_NoError;
    errcode = readTLE7209(chipSelectPin, TLE7209_READ_IDENTIFIER, &deviceID);

    if(errcode){
        #if DEBUG
            Serial.println("TLE7209readIDandVersion() failed");
        #endif
        return errcode;
    }

    // Chip ID is fixed, check if match
    if(deviceID != TLE7209_DEFAULT_DEVICE_ID){
        #if DEBUG
            Serial.println("Wrong Device ID detected\n");
        #endif
        return TLE7209_InvalidDeviceId;
    }

    errcode = readTLE7209(chipSelectPin, TLE7209_READ_VERSION, version);
    return errcode;
}


TLE7209_Error TLE7209readDiagnosticRegister(byte chipSelectPin,
                                            byte *diagnostics){                 // TODO: docstring

    TLE7209_Error errcode = TLE7209_NoError;
    errcode = readTLE7209(chipSelectPin,
                          TLE7209_READ_DIAG_REGISTER,
                          diagnostics);

    if(errcode){
        #if DEBUG
            Serial.println("TLE7209readDiagnosticRegister() failed");
        #endif
        return errcode;
    }

    // All bits except the MSB are set to 1 when no error occurred;
    // The MSB is just for info: mask out the MSB, then check against
    // the no-error condition
    if ((*diagnostics & TLE7209_ALL_ERROR_BITS) != TLE7209_ALL_ERROR_BITS){
        #if DEBUG
            Serial.println("DIA_REG: several bits set\n");
        #endif
        return TLE7209_DiagnosticsError;
    }
    return TLE7209_NoError;
}


// TODO: may implement an extra function that gives more info about
// which diagnostics error happened, e.g., readErrorDetails
