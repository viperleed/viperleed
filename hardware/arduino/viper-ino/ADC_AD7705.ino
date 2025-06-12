/*
ViPErLEED - Driver for AD7705 Analog-to-Digital Converter
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 26.04.2021
---------------------
*/


#include "ADC_AD7705.h"


/* ---------- AD7705 ADC FUNCTIONS ---------- */

/** Starts an I/O operation for the AD7705 with given chip select pin */
void AD7705startIO(byte chipSelectPin) {
    SPI.beginTransaction(AD7705_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
}

/** Finishes an I/O operation for the AD7705 with given chip select pin */
void AD7705endIO(byte chipSelectPin) {
    digitalWrite(chipSelectPin, HIGH);
}

/** Resets I/O of the AD7705 by writing 32 high bits */
void AD7705resetCommunication(byte chipSelectPin) {
    AD7705startIO(chipSelectPin);
    SPI.transfer16(0xffff);
    SPI.transfer16(0xffff);
    AD7705endIO(chipSelectPin);
}

/** Writes the clock register of the AD7705 with the given chip select pin
 *  The update rate can be AD7705_50HZ, AD7705_60HZ, or AD7705_500HZ
 *  Not used because setting the update rate requires self-calibration */
void AD7705setClock(byte chipSelectPin, byte updateRate) {
   AD7705startIO(chipSelectPin);
   SPI.transfer(AD7705_REG_CLOCK);
   SPI.transfer(AD7705_CLK | updateRate);
   AD7705endIO(chipSelectPin);
 }

void AD7705setGainAndTrigger(byte chipSelectPin, byte channel, byte gain) {
    /**Set the gain, and trigger an AD7705 ADC.

    This essentially sets a well-defined point in time after
    which conversions of the selected input channel will be
    available, after amplification with the selected gain.

    Notice that this function does not initiate the reading
    of any of the conversion results from the Arduino. Thus,
    every 1/updateRate sec after this, new values will be
    available in the data register, and they will be thrown
    away if they are not read.

    The function also sets bipolar unbuffered operation mode,
    as well as 'normal operation' (i.e., not a calibration)
    that is necessary for our inputs. Any old data that may
    still be present in the data register is also discarded.

    Parameters
    ----------
    chipSelectPin : byte
        The Arduino pin that corresponds to the ADC that has to
        be addressed
    channel : byte, {0, 1}
        Which of the input channels will be used for conversion
    gain : byte, {0...7}
        Which gain value should be used. The input is amplified
        by a factor 2**gain.
    **/
    AD7705startIO(chipSelectPin);
    //SPI.transfer16(0xffff); //reset should not be necessary
    //SPI.transfer16(0xffff);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_TRIGGER | gain << 3);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(gain << 3);

    // Read the data register to discard any old data
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0xffff);
    AD7705endIO(chipSelectPin);
}

void AD7705selfCalibrate(byte chipSelectPin, byte channel,
                         byte gain, byte updateRate) {
    /**Self-calibrate an AD7705 ADC.

    Thereafter, one has to wait until self-calibration is
    done, by calling AD7705waitForCalibration, or by waiting
    for the first data (AD7705waitAndReadData).

    Triggering during self-calibration should never be performed.

    Parameters
    ----------
    chipSelectPin : byte
        The Arduino pin corresponding to the ADC
        that needs to be calibrated.
    channel : byte, {0, 1}
        Which of the input channels needs to be
        calibrated.  Only this channel will be
        calibrated. The calibration of the other
        channel remains untouched.
    gain : byte, {0...7}
        Which gain should be used for calibration.
        The signal will be amplified by a factor
        2**gain before conversion.  A gain switch
        always requires a new calibration.
    updateRate : {AD7705_50HZ, AD7705_60HZ, AD7705_500HZ}
        The input conversion rate for which calibration
        needs to be done. When changing updateRate, the
        ADC always needs a new calibration.
    **/
    AD7705startIO(chipSelectPin);
    SPI.transfer(AD7705_REG_CLOCK | channel);
    SPI.transfer(AD7705_CLK | updateRate);     // Set clock and update rate
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_SELFCAL | gain << 3);  // Start self-calibration
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0xffff);                    // Read data register to ensure DRDY is off.
    AD7705endIO(chipSelectPin);
}

int32_t AD7705getCalibrationRegister(byte chipSelectPin, byte channel,
                                     byte calibrationRegister) {                // TODO: Datasheet says one should write high to FSYNC before doing operations with the calibration registers
    /**Return the offset or gain calibration register of an AD7705.

    This function can be called after a call to
    AD7705selfCalibrate() to retrieve the calibration
    result, with a short delay(1) in between.

    Parameters
    ----------
    chipSelectPin : byte
        The Arduino pin acting as chip-select for the
        ADC from which the register has to be read
    channel : byte
        The ADC channel whose calibration register
        has to be read
    calibrationRegister : {AD7705_REG_OFFSET, AD7705_REG_GAIN}
        The register whose contents are to be read.
        No check is done that in fact calibrationRegister
        is one of the admissible values!                                        // TODO: perhaps we should check this
    */
    AD7705startIO(chipSelectPin);
    SPI.transfer(calibrationRegister | AD7705_READ_REG | channel);
    int32_t result = (int32_t)SPI.transfer16(0xffff) << 8; //read high 16 bits
    result |= SPI.transfer(0xff);              //read low 8 bits
    AD7705endIO(chipSelectPin);
    return result;
}

// TODO: why is this returning stuff? Should be void.
/** Sets the value in an AD7705 offset or gain calibration register.
 *  'theRegister' must be AD7705_REG_OFFSET or AD7705_REG_GAIN.
 *  This function can be called to restore the result of a previous self-calibration */
int32_t AD7705setCalibrationRegister(byte chipSelectPin, byte channel, byte theRegister, int32_t value) {
  AD7705startIO(chipSelectPin);                                                 // TODO: Datasheet says one should write high to FSYNC before doing operations with the calibration registers
  SPI.transfer(theRegister | channel);
  SPI.transfer16(value>>8);                  //write high 16 bits
  SPI.transfer(value&0xff);                  //write low 8 bits
  AD7705endIO(chipSelectPin);
}

/** Waits until the AD7705 with the given chip select pin has finished self-calibration */
void AD7705waitForCalibration(byte chipSelectPin, byte channel) {
    uint32_t timeoutCounter = 0;
    AD7705startIO(chipSelectPin);
    do {
        if (timeoutCounter > 40100) //40100 == 4.01 seconds > than our arduino timeout (4 seconds)
            return;
        delayMicroseconds(AD7705_DELAY_MICRO);
        SPI.transfer(AD7705_REG_SETUP | AD7705_READ_REG | channel);
        timeoutCounter++;
    } while (SPI.transfer(0x0) & AD7705_SELFCAL);  //read setup register and check for self-calibration
    AD7705endIO(chipSelectPin);
}

/** Waits for data and then reads the given channel of the AD7705 with the given chip select pin.
 *  Returns a signed 16-bit int with the signed value (N.B. we use bipolar mode only) */
int16_t AD7705waitAndReadData(byte chipSelectPin, byte channel) {
    uint32_t timeoutCounter = 0;
    while (AD7705readCommRegister(chipSelectPin, channel) & AD7705_DRDY) {//DRDY bit is 0 when ready
        if (timeoutCounter > 40100) //40100 == 4.01 seconds > than our arduino timeout (4 seconds)
            return -1;
        delayMicroseconds(AD7705_DELAY_MICRO);
        timeoutCounter++;
    }
    AD7705startIO(chipSelectPin);
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    int result = SPI.transfer16(0xffff);
    result ^= 0x8000; //flip highest bit: AD7705 output is 0 for -Vref, 0x8000 for 0, 0xffff for +Vref
    AD7705endIO(chipSelectPin);
    return result;
}

/** Reads and returns the communications register of the AD7705 with given chip select pin.
 *  This operation should be done for the channel currently in use. */
byte AD7705readCommRegister(byte chipSelectPin, byte channel) {
    AD7705startIO(chipSelectPin);
    //SPI.transfer16(0xffff); //reset should not be necessary
    //SPI.transfer16(0xffff);
    SPI.transfer(AD7705_READ_REG | AD7705_REG_COMM | channel);
    byte result = SPI.transfer(0xff);
    AD7705endIO(chipSelectPin);
    return result;
}
