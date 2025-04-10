/* ViPErLEED - Firmware for Arduino hardware controller

@author: Bernhard Mayr
@author: Michael Schmid (@schmid-iap)
@author: Michele Riva (@michele-riva)
@author: Florian Dörr (@FlorianDoerr)
Date: 26.04.2021
*/

#ifndef _VIPERLEED
#define _VIPERLEED


// arduino_utils.h, states-def.h and viper-serial.h come from ../lib
#include "arduino_utils.h"  // from ../lib; for setChipSelectHigh,
                            // getMedian, bigger, biggest
#include "states-def.h"     // Basic state-manchine definitions
#include "viper-serial.h"   // Serial-communication functions and constants
#include "ADC_AD7705.h"     // Settings of the currently used ADC
#include "DAC_AD5683.h"     // Settings of the currently used DAC


// Useful macros, structs and unions
#define LENGTH(array)  (sizeof(array) / sizeof((array)[0]))
#define MIN(a, b)      ((a) < (b) ? (a) : (b))
#ifndef NAN
    #define NAN   (1.0/0)   // Floating-point not-a-number
#endif

// Interpret the same value as either uint16_t or a 2-bytes array
union uint16OrBytes{
    uint16_t asInt;
    byte asBytes[2];
};

// Interpret the same value as either int16_t or a 2-bytes array
union int16OrBytes{
    int16_t asInt;
    byte asBytes[2];
};

// Interpret the same value as either int32_t or a 4-bytes array
union int32OrBytes{
    int32_t asInt;
    byte asBytes[4];
};

// Interpret the same value as either (32-bit) float or a 4-bytes array
union floatOrBytes{
    float asFloat;
    byte asBytes[4];
};


/** ------------------------ Communication with PC ----------------------- **/
// Further globals and constans defined in viper-serial.h

// Acceptable messages for communication with the PC

// PC_AUTOGAIN: PC requested auto-gain for ADCs (ASCII 'A').
#define PC_AUTOGAIN          65
// PC_CALIBRATION: PC requested self-calibration
// of all ADCs at all gains (ASCII 'C').
#define PC_CALIBRATION       67
// PC_CHANGE_MEAS_MODE: PC requested a change between
// continuous and single-measurement modes (ASCII 'm').
#define PC_CHANGE_MEAS_MODE 109
// PC_MEASURE_ONLY: PC requested measurement
// without changing voltage (ASCII 'M').
#define PC_MEASURE_ONLY      77
// PC_SET_UP_ADCS: PC requested to prepare the
// ADCs for a measurement (ASCII 'S').
#define PC_SET_UP_ADCS       83                                                 // TODO: The python side will have to keep track of when the last calibration was done, and warn if it is too old.
// PC_SET_VOLTAGE: PC requested to set a certain energy (ASCII 'V').
// After setting the voltage, the arduino will start measuring.
#define PC_SET_VOLTAGE       86
// PC_SET_VOLTAGE_ONLY: PC requested set energy without
// follow-up measurement (ASCII 'v').
#define PC_SET_VOLTAGE_ONLY 118

// Serial-communication-related error codes
// Further errors defined in viper-serial.h and states-def.h

// ERROR_NEVER_CALIBRATED: The ADCs have never
// been calibrated before since bootup.
#define ERROR_NEVER_CALIBRATED    6
// ERROR_ADC_SATURATED: One of the ADC values reached saturation,
// and gain can't be decreased further.
#define ERROR_ADC_SATURATED       8
// ERROR_TOO_HOT: The temperature read by the LM35 is too high.
#define ERROR_TOO_HOT             9
// ERROR_HARDWARE_UNKNOWN: The PC never asked for the hardware configuration.
#define ERROR_HARDWARE_UNKNOWN   10

// Further globals and constans defined in states-def.h

// STATE_SET_UP_ADCS: Pick correct ADC channels and no. of measurement points.
#define STATE_SET_UP_ADCS          1
// STATE_SET_VOLTAGE: Set a voltage with the DAC, wait, then trigger the ADCs.
#define STATE_SET_VOLTAGE          2
// STATE_CHANGE_MEASUREMENT_MODE: Get and set the desired measurement mode.
#define STATE_CHANGE_MEASUREMENT_MODE 3
// STATE_MEASURE_ADCS: ADC measurements in progress.
#define STATE_MEASURE_ADCS         4
// STATE_ADC_VALUES_READY: ADC measurements done.
#define STATE_ADC_VALUES_READY     5
// STATE_AUTOGAIN_ADCS: Find optimal gain for both ADCs.
#define STATE_AUTOGAIN_ADCS        6   7
// STATE_CALIBRATE_ADCS: Figure out correct offset and
// calibration factors for ADCs at all gains.
#define STATE_CALIBRATE_ADCS       8
// continuousMeasurement: Decides if the Arduino continues to measure
// and return data or if it stops after doing so once.
bool continuousMeasurement = false;


/** --------------------- Hardware-specific settings --------------------- **/
// Measurement devices

// N_MAX_MEAS: Number of measurement devices that we can have: 3
// ADC#0, ADC#1, LM35 (if present)
#define N_MAX_MEAS         3
// N_MAX_ADCS_ON_PCB: Maximum number of ADCs on the PCB
#define N_MAX_ADCS_ON_PCB  2

// I/O pins on Arduino board

// CS_DAC: (Inverted) chip select of DAC.
#define CS_DAC             4
// CS_ADC_0: (Inverted) chip select of ADC #0.
#define CS_ADC_0           5
// CS_ADC_1: (Inverted) chip select of ADC #1 (optional).
#define CS_ADC_1           6
// LM35_PIN: LM35 temperature sensor (optional).
#define LM35_PIN          A1
// RELAY_PIN: Relay is connected here if present.
#define RELAY_PIN         A4
// JP_I0_PIN: Jumper for manual I0 2.5V range, same as relay pin.
#define JP_I0_PIN         A4
// JP_AUX_PIN: Jumper for manual AUX 2.5 V range.
#define JP_AUX_PIN        A5

// Arduino internal reference, ADC maximum, and settings for relay and LM35

// VREF_INTERNAL: Arduino micro internal ADC reference is 2.56 V.
#define VREF_INTERNAL   2.56
// ARDUINO_ADC_MAX: 10-bit internal ADC
#define ARDUINO_ADC_MAX  0x3ff
// RELAY_MIN_ADU: Relay has about 0.12 V with pull-up,
// ADC signal is larger than this.
#define RELAY_MIN_ADU     24
// RELAY_MAX_ADU: ADC signal of relay with pull-up is less than this.
#define RELAY_MAX_ADU     96
// LM35_MAX_ADU: LM35 signal should be less than 0.8 V (80 degC).
#define LM35_MAX_ADU     320
// ARDUINO_ADU_VOLTS: Internal ADC: volts per ADU.
#define ARDUINO_ADU_VOLTS (VREF_INTERNAL/ARDUINO_ADC_MAX)
// ARDUINO_ADU_DEG_C: LM35 degC per ADU.
#define ARDUINO_ADU_DEG_C (100*ARDUINO_ADU_VOLTS)

// Scaling of ADC input ranges due to voltage
// dividers or preamplification on the board.

// ADC_0_CH0_SCALE_JO: ADC#0 channel 0:
// with JP at I0 open or relay off, in volts
#define ADC_0_CH0_SCALE_JO      0.004
// ADC_0_CH1_SCALE: ADC#0 channel 1: High-voltage divider
#define ADC_0_CH1_SCALE    (16000./39.*0.001)
// ADC_1_CH0_SCALE: ADC#1 channel 0: uAmps
#define ADC_1_CH0_SCALE        (10/2500.0)
// ADC_1_CH1_SCALE_JO: ADC#1 channel 1:
// with JP5 open (voltage divider), in volts
#define ADC_1_CH1_SCALE_JO      0.004

// Which hardware and closed jumpers do we have?
// Bits of present/closed jumpers are 1.

// ADC_0_PRESENT: Bit is set if ADC #0 was detected
#define ADC_0_PRESENT   0x01
// ADC_1_PRESENT: Bit is set if ADC #1 was detected
#define ADC_1_PRESENT   0x02
// LM35_PRESENT: Bit is set if LM35 temperature sensor was detected
#define LM35_PRESENT    0x04
// RELAY_PRESENT: Bit is set if the relay for I0 input 2.5V/10V is present
#define RELAY_PRESENT   0x08
// JP_I0_CLOSED: Bit is set if JP3 7-8 is closed
// or relay on (to indicate 2.5V I0 range)
#define JP_I0_CLOSED    0x10
// JP_AUX_CLOSED: Bit is set if JP5 is closed
// (to indicate 2.5V AUX range rather than 10V range)
#define JP_AUX_CLOSED   0x20

// ADC saturation thresholds

// ADC_SATURATION: Threshold to deem abs(adcValue) at saturation,
// in principle it should be 32767/-32768 but needs testing.
#define ADC_SATURATION       32760
// ADC_RANGE_THRESHOLD: If output is larger than this
// (~25% of range), we need a new gain next time
#define ADC_RANGE_THRESHOLD  0x3fff


/** ------------------- Globals for firmware functions ------------------- **/
// Timers (defined in milliseconds)

// dacSettlingTime: The time interval for the DAC output to be stable.
// This is just a default value. The actual one comes from the PC with a
// PC_SET_VOLTAGE command, and is read in setVoltage().
uint16_t      dacSettlingTime = 100;
// nextVoltageStep: Counter for multiple voltage steps.
byte nextVoltageStep;

// hardwareDetected: Bits set indicate this hardware is present/jumper closed.
uint16OrBytes hardwareDetected;
// hardwareNeverChecked: Is set false if the
// PC asked for the hardware configuration.
bool hardwareNeverChecked = true;
// takeMeasurements: Is set according to set voltage command.
// True if measurements are requested.
bool takeMeasurements = true;

// ADCs: measurement frequency, channel, gain

// adcUpdateRate: Update rate for both ADCs (will be set for line frequency).
byte     adcUpdateRate;
// adc0Channel: Channel of ADC#0
byte     adc0Channel;
// adc1Channel: Channel of ADC#1
byte     adc1Channel;
// adc0Gain: Gain of ADC#0, 0...7
byte     adc0Gain = 0;
// adc1Gain: Gain of ADC#1, 0...7
byte     adc1Gain = 0;
// adc0ShouldDecreaseGain: Whether ADC#0 should
// increase its gain in the next cycle.
bool     adc0ShouldDecreaseGain = false;
// adc1ShouldDecreaseGain: Whether ADC#1 should
// increase its gain in the next cycle.
bool     adc1ShouldDecreaseGain = false;

// ADCs: quantities needed for self-calibration

// calibrationGain: Gain for which a calibration is currently
// being performed (in parallel for both ADCs).
byte     calibrationGain = 0;
// calibratedChannels: Array of flags to keep track of which
// channels (last index) of the ADCs have been calibrated.
bool     calibratedChannels[N_MAX_ADCS_ON_PCB][2];
// selfCalDataVsGain: For each gain, two ADCs, two
// channels each, and last index is offset(0)&gain(1).
int32_t  selfCalDataVsGain[AD7705_MAX_GAIN + 1][2][2][2];

// ADCs: variables for measuring and storing the line frequency ripple

// maximumPeak: Maximum of measurement, one for each ADC.
int16_t  maximumPeak[N_MAX_ADCS_ON_PCB];
// : minimumPeak: Minimum of measurement, one for each ADC.
int16_t  minimumPeak[N_MAX_ADCS_ON_PCB];
// adc0RipplePP: Ripple (peak-peak) measured
// at gain=0 for ADC#0 during auto-gain.
int16_t  adc0RipplePP = 0;
// adc1RipplePP: Ripple (peak-peak) measured
// at gain=0 for ADC#1 during auto-gain.
int16_t  adc1RipplePP = 0;

/* // ADC container                                                             // TODO: use it to simplify code below
struct analogToDigitalConverter {
    byte     chipSelect;    // to be initialized!
    uint16_t present;       // to be initialized!
    byte     channel = AD7705_CH0;
    byte     gain = 0;
    bool     shouldDecreaseGain = false;
    // Ripple (peak-peak) measured at gain=0 during auto-gain.
    int16_t  ripplePP = 0;                                                       // TODO: should we keep the largest ripple value ever measured or the latest?
    // For each gain, two channels each, and last index is offset(0) & gain(1)
    int32_t  calibrationForGain[AD7705_MAX_GAIN + 1][2][2];
    bool     calibratedChannels[2] = {false, false};
} externalADCs[N_MAX_ADCS_ON_PCB]; */

// Measurements

// numMeasurementsToDo: No. of ADC measurements
// to do before measurement is considered over.
uint16_t numMeasurementsToDo = 1;
// numMeasurementsToDoBackup: Copy of the previous one, used
// to restore the previous value after auto-gain is done.
uint16_t numMeasurementsToDoBackup = 1;
// numMeasurementsDone: Counter for the number
// of ADC measurements done so far.
uint16_t numMeasurementsDone;
// summedMeasurements: Measurements of ADCs and LM35 are summed up here.
int32_t  summedMeasurements[N_MAX_MEAS];
// continuousMeasurementInterval: Time between measurements
// if continuous-measurement mode is on. Currently unused.
uint16_t continuousMeasurementInterval;

// floatOrBytes: Measurements in physical units
floatOrBytes fDataOutput[N_MAX_MEAS];                                           // TODO: rename

#endif
