/*
ViPErLEED - Driver for TLE7209-3R H-Bridge
---------------------
Authors: Michele Riva, Christoph Pfungen
21.04.2023, MR, CP: First version
26.11.2024, SM: Overhaul whole code as class; only use SPI, when needed
---------------------
*/

#ifndef _VIPERLEED_TLE7209
#define _VIPERLEED_TLE7209

#define TLE7209_DELAY_MICRO 100

#define TLE7209_USE_SPI     false

#if TLE7209_USE_SPI
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


// TLE7209 verification byte definitions
#define TLE7209_TRANS_F 0x01	// Bit is set if previous transfer was recognized as valid
                                // Bit is cleared on error during previous transfer

#endif

// TLE7209 error codes
enum TLE7209_Error {
    TLE7209_NoError,
    TLE7209_TransmissionError,
    TLE7209_InvalidDeviceId,
    TLE7209_DiagnosticsError,
    TLE7209_SPINotAvailable,
};


class TLE7209 : public MotorDriver {
    public:
        TLE7209(byte enable_pin)
        : MotorDriver(enable_pin) {}

        TLE7209(byte chip_select_pin, byte enable_pin)
        : MotorDriver(chip_select_pin, enable_pin) {}

        void setup() override {
#if TLE7209_USE_SPI
            if (_SPI_CS_PIN) {
                setChipSelectHigh(_SPI_CS_PIN);
            }
#endif
            pinMode(_ENABLE_PIN, OUTPUT);    // Use for TLE7209 enable pin 'EN'
            digitalWrite(_ENABLE_PIN, HIGH);  // 'EN' == 1: Activate TLE7209 output drivers
        }

        TLE7209_Error get_version(byte* version) {
            return TLE7209readIDandVersion(_SPI_CS_PIN, version);
        }

        TLE7209_Error get_diagnostic_info(byte* info) {
            return TLE7209readDiagnosticRegister(_SPI_CS_PIN, info);
        }

    private:
#if TLE7209_USE_SPI
        /** Start an I/O operation for the TLE7209 with given chip-select pin.**/
        void TLE7209startIO(byte chipSelectPin) {
            SPI.beginTransaction(TLE7209_SPI_SETTINGS);
            digitalWrite(chipSelectPin, LOW);
        }

        /** Finish an I/O operation for the TLE7209 with given chip-select pin.**/
        void TLE7209endIO(byte chipSelectPin) {
            digitalWrite(chipSelectPin, HIGH);
            SPI.endTransaction();
        }
#endif

        TLE7209_Error readTLE7209(byte chipSelectPin, byte request, byte *data) {
            /**Generic SPI routine to read a register from the TLE7209.

            Parameters
            ----------
            chipSelectPin : byte
                Pin on the Arduino Micro acting as the chip select line for the TLE7209
            request : {TLE7209_READ_IDENTIFIER, TLE7209_READ_VERSION,
                      TLE7209_READ_DIAG_REGISTER}
                Describes the register or value which should be retrieved
            data : byte *
                After a successful read, this is where the register content is stored

            Returns
            -------
            TLE7209_Error : enum
                TLE7209_NoError for successful read
                TLE7209_TransmissionError if verification byte indicates errors
            **/
#if TLE7209_USE_SPI
            if (!_SPI_CS_PIN) {
                return TLE7209_SPINotAvailable;
            }
            TLE7209startIO(chipSelectPin);
            uint16_t bytesRead = SPI.transfer16((uint16_t)(request << 8));
            TLE7209endIO(chipSelectPin);

            // Check verification byte
            uint8_t transmit_ok = bytesRead >> 8;
            transmit_ok &= 0b00111111; // The highest two bits are not relevant
            if(transmit_ok != TLE7209_SPI_TRANSMISSION_OK) {
                #if DEBUG
                    Serial.println("Verification byte: TRANS_F is set and/or wrong bit toggle sequence detected\n");
                #endif
                return TLE7209_TransmissionError;
            }

            // Return only the LSB (data byte)
            *data = (uint8_t)bytesRead;
            return TLE7209_NoError;
#else
            return TLE7209_SPINotAvailable;
#endif
        }

        // Perform two transfers back-to-back (READ_ID, READ_VER)
        TLE7209_Error TLE7209readIDandVersion(byte chipSelectPin, byte *version) {
            /**Read device ID and chip revision from the TLE7209.

            Parameters
            ----------
            chipSelectPin : byte
                Pin on the Arduino Micro acting as the chip select line for the TLE7209
            version : byte *
                After a successful read, this is where the chip revision is stored

            Returns
            -------
            TLE7209_Error : enum
                TLE7209_NoError for successful read
                TLE7209_TransmissionError if verification byte indicates errors
                TLE7209_InvalidDeviceId if device ID does not match fixed value,
                cf. datasheet sec. 2.4.2.7 (p. 17)
            **/
#if TLE7209_USE_SPI
            byte deviceID;
            TLE7209_Error errcode = TLE7209_NoError;
            errcode = readTLE7209(chipSelectPin, TLE7209_READ_IDENTIFIER, &deviceID);

            if(errcode) {
                #if DEBUG
                    Serial.println("TLE7209readIDandVersion() failed");
                #endif
                return errcode;
            }

            // Chip ID is fixed, check if match
            if(deviceID != TLE7209_DEFAULT_DEVICE_ID) {
                #if DEBUG
                    Serial.println("Wrong Device ID detected\n");
                #endif
                return TLE7209_InvalidDeviceId;
            }

            errcode = readTLE7209(chipSelectPin, TLE7209_READ_VERSION, version);
            return errcode;
#else
            return TLE7209_SPINotAvailable;
#endif
        }

        TLE7209_Error TLE7209readDiagnosticRegister(byte chipSelectPin,
                                                    byte *diagnostics) {
            /**Read diagnostics from the TLE7209.
            For further details see datasheet section 2.4.2.6.

            Parameters
            ----------
            chipSelectPin : byte
                Pin on the Arduino Micro acting as the chip select line for the TLE7209
            diagnostics : byte *
                After a successful read, this is where the register content is stored

            Returns
            -------
            TLE7209_Error : enum
                TLE7209_NoError for successful read
                TLE7209_TransmissionError if verification byte indicates errors
                TLE7209_DiagnosticsError if one or more status flags have been
                set by the TLE7209

            **/
#if TLE7209_USE_SPI
            TLE7209_Error errcode = TLE7209_NoError;
            errcode = readTLE7209(chipSelectPin,
                                  TLE7209_READ_DIAG_REGISTER,
                                  diagnostics);
            if(errcode) {
                #if DEBUG
                    Serial.println("TLE7209readDiagnosticRegister() failed");
                #endif
                return errcode;
            }

            // All bits except the MSB are set to 1 when no error occurred;
            // The MSB is just for info: mask out the MSB, then check against
            // the no-error condition
            if ((*diagnostics & TLE7209_ALL_ERROR_BITS) != TLE7209_ALL_ERROR_BITS) {
                #if DEBUG
                    Serial.println("DIA_REG: several bits set\n");
                #endif
                return TLE7209_DiagnosticsError;
            }
            return TLE7209_NoError;
#else
            return TLE7209_SPINotAvailable;
#endif
        }

};


#endif  // _VIPERLEED_TLE7209
