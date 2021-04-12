/*
Main LeedControl File
---------------------
Author: Bernhard Mayr, Michael Schmid, Michele Riva, Florian Dörr
Date: 11.04.2021
---------------------
*/
//Libraries
#include <Arduino.h>
#include <SPI.h>

//File with Arduino settings
//#include "viper-ino.h"
//#include "viper-ino_control.h"

#define DEBUG           true  //debug mode, writes to serial line, for use in serial monitor

//Useful macros
#define LENGTH(array)  (sizeof(array) / sizeof((array)[0]))
#ifndef NAN                                                                                            // !!!!!!!!!! This is never used
  #define NAN        (1.0/0)   //floating-point not-a-number
#endif

//global variables for communication with Python module
#define STARTMARKER 254
#define ENDMARKER 255
#define SPECIAL_BYTE 252                                         // !!!!!!!!!!!!!!!!!!!!!! CHANGE FOR PYTHON (Python has value 253)
#define PC_ERROR 253                                             // !!!!!!!!!!!!!!!!!!!!!! Define in Python
#define MAX_MSG_LENGTH 16

//definition of messages for communication with PC
#define PC_HARDWARE 3
#define PC_INIT_ADC 4
#define PC_OK 5
#define PC_MEASURE 6
#define PC_SET_VOLTAGE 7
#define PC_AUTOGAIN 8
#define PC_RESET 82                                               // !!!!!!!!!!!!!!!!!!!!!! CHANGE FOR PYTHON (Python has value 68)
#define PC_DEBUG 0

// << THINGS TO DO EVERYWHERE (to conform with c++ style):
// - use "lower-camel-case" for variables (i.e., thisIsANewVariable) instead of "lower-snake-case" (i.e., this_is_a_new_variable)
// - use "upper-camel-case" for constants (i.e., THIS_IS_A_CONSTANT). I think most are already correct
// - use "lower-camel-case" for functions

//Variables for the communication protocoll
byte numBytesRead = 0;            //counter for received bytes
byte msgLength = 0;           //numbers of bytes received in message, 2nd byte of received message !!!!!!!!!!!!!!!!!!!!!!!!!!!!   This is set in readFromSerial but never used outside of it. Perhaps we should use it to check that indeed we got the data that we expected.
byte data_received_counter = 0;     //total number of bytes received  !!!!!!!!!!!!!!!!!!  This is also used only locally in decodeMessage, and it's probably better to keep it only local. There should be minimal to no performance hit by doing this. Also, renaming it to something clearer would help
byte data_received[MAX_MSG_LENGTH];    //Real message  !!!!!!!!!!! Judging from the code, and if I didn't misinterpret anything, I think we can have only one buffer for the input/output to the PC, so we can replace "data_received" and "data_send" with a single "message" array
byte data_send[MAX_MSG_LENGTH];  
byte inputRegister[MAX_MSG_LENGTH];

//Finite state machine 
int currentState = 0; //!!!!!!!!!!!!!  Also, perhaps it would be better to have them #define(d) rather than as const int? (save a little memory, 2 bytes * 7 names)
const int STATE_IDLE = 0; 
const int STATE_SETUP_ADC = 1;
const int STATE_SET_VOLTAGE = 2;
const int STATE_TRIGGER_ADCS = 3; 
const int STATE_ADC_MEASURE = 4; 
const int STATE_ADC_VALUE_READY = 5;
const int STATE_AUTOGAIN = 6; 
const int STATE_GET_HARDWARE = 7; 
const int STATE_ERROR = 8; 

//Flags 
boolean readingFromSerial = false;            //Flag True when arduino receives new data
boolean newMessage = false;               //Fag True when all new data has arrived

//ADC settings
#define AD7705_MAX_GAIN 7

//ADC measurement settings
uint16_t numMeasurementsToDo = 1;
uint16_t numMeasurementsDone = 0;

//ADC saturation markers
#define ADC_POSITIVE_SATURATION 0xffff
#define ADC_NEGATIVE_SATURATION 0x0000
#define ADC_RANGE_THRESHOLD 0x3fff

//Timers (defined in milliseconds)
unsigned long initialTime; 
unsigned long timeout = 5000;      // !!!!!!!!!!!! this can become a #define, or at least a const    
unsigned long dac_settletime = 100; // !!!!!!!!!! not sure why we use a long here, as we accept only 2 bytes (uint16_t) in setVoltage()
unsigned long currentTime;

//Number of measurement devices that we can have: ADC#0, ADC#1, LM35 (if present)
#define N_MAX_MEAS         3
//Arduino I/O pins
#define CS_DAC             4   //(inverted) chip select of DAC
#define CS_ADC_0           5   //(inverted) chip select of ADC #0
#define CS_ADC_1           6   //(inverted) chip select of ADC #1 (optional)
#define LM35_PIN          A1   //LM35 temperature sensor (optional)
#define RELAY_PIN         A4   //relay is connected here if present
#define JP_I0_PIN         A4   //jumper for manual I0 2.5V range, same as relay pin
#define JP_AUX_PIN        A5   //jumper for manual AUX 2.5 V range
#define WORKINPROGRESS     7   //Python script will decide which pin !!!!!!!!!!! this needs to be NOT a #define because it needs updating. It should be a byte. Also, probably good to give it a more descriptive name, like, whatToMeasure. with 8 bits, we can have a bit representation of the devices. E.g.: ADC0 if WORKINPROGRESS == 0x01 == 0b00000001, ADC1 if WORKINPROGRESS == 0x02 == 0b00000010, LM35 if WORKINPROGRESS == 0x04 == 0b00000100. And combinations can be obtained by bitwise OR. In principle we could also encode it better, such that this single byte can be used also to determine which channel of which ADC to measure: e.g.: MEASURE_ADC0_CH0 = 0x01, MEASURE_ADC0_CH1 = 0x02, MEASURE_ADC1_CH0 = 0x04, MEASURE_ADC1_CH1 = 0x08, MEASURE_LM35 = 0x0f. Obviously then one needs to handle somehow the cases in which the python asks for measuring two channels of the same ADC, which cannot be done at the same time.

//Arduino internal reference voltage, ADC maximum, etc
#define VREF_INTERNAL    2.56  //arduino micro internal ADC reference is 2.56 V
#define ARDUINO_ADC_MAX  0x3ff //10-bit ADC
#define ARDUINO_ADU_VOLTS (VREF_INTERNAL/ARDUINO_ADC_MAX)      //volts per ADU
#define ARDUINO_ADU_DEG_C (100*VREF_INTERNAL/ARDUINO_ADC_MAX)  //LM35 degC per ADU
#define RELAY_MIN_ADU     24   //relay has about 0.12V with pullup, ADC signal is larger than this
#define RELAY_MAX_ADU     96   //ADC signal of relay with pullup is less than this
#define LM35_MAX_ADU     320   //LM35 signal should be less than 0.8 V (80 degC)

// Definitions for AD7705 analog-to-digital converter
#define AD7705_SPIMODE     3   //data accepted at rising edge of SCLK
#define AD7705_MAX_GAIN    7   //gains are 0 (x1) to 7 (x128)
//AD7705 communication register
#define AD7705_READ_REG  0x08  //set bit for reading a register
#define AD7705_CH0       0x00  //set bit for adressing channel 0
#define AD7705_CH1       0x01  //set bit for adressing channel 1
#define AD7705_REG_COMM  0x00  //bits to set for accessing communication register
#define AD7705_REG_SETUP 0x10  //bits to set for accessing setup register
#define AD7705_REG_CLOCK 0x20  //bits to set for accessing clock register
#define AD7705_REG_DATA  0x30  //bits to set for accessing data register (16 bit)
#define AD7705_REG_OFFSET 0x60 //bits to set for accessing calibration offset register (24 bit)
#define AD7705_REG_GAIN  0x70  //bits to set for accessing calibration gain register (24 bit)
#define AD7705_DRDY      0x80  //bit mask for Data Ready bit in communication register
//AD7705 setup register
#define AD7705_STATE_TRIGGER_ADCS     0x01  //bit to set for FSYNC in setup register
#define AD7705_BUFMODE   0x02  //bit to set for buffered mode (not used by us)
#define AD7705_UNIPOLAR  0x04  //bit to set for unipolar mode (not used by us)
#define AD7705_SELFCAL   0x40  //bit to set for self-calibration in setup register
//AD7705 clock register
#define AD7705_CLK       0x0C  //bits to set for 4.9152 MHz crystal: CLK=1, CLKDIV=1
#define AD7705_50HZ      0x04  //bits to set for 50 Hz update rate and suppression
#define AD7705_60HZ      0x05  //bits to set for 60 Hz update rate and suppression
#define AD7705_500HZ     0x07  //bits to set for 500 Hz update rate
// Correction for the finite source impedance and ADC input impedance
#define R_SOURCE       1300.0  //source resistance of our circuit at ADC inputs, 1.3 kOhm
//#define R_IN_0     (1.0/(7e-12*38400)) //input resistance of ADC7705 in gain0 = x1, nominal (3.72MegOhm)
#define R_IN_0          4.1e6  //input resistance of ADC7705 in gain0 = x1, measured typical value
#define REF_OVER_RANGE  (2500.0/32768)    //millivolts (uA@1kOhm) per bit at gain0, bipolar, input 0ohm
const float voltsPerBit[] = {
        REF_OVER_RANGE/R_IN_0*(R_IN_0 + R_SOURCE),         //gain0 = x1
        REF_OVER_RANGE/R_IN_0*2*(R_IN_0/2 + R_SOURCE)/2,   //gain1 = x2 has half R_in_0
        REF_OVER_RANGE/R_IN_0*4*(R_IN_0/4 + R_SOURCE)/4,   //gain2 = x4 has 1/4 R_in_0
        REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/8,   //gain3 = x8 and up: 1/8 R_in_0
        REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/16,  //gain4 = x16
        REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/32,  //gain5 = x32
        REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/64,  //gain6 = x64
        REF_OVER_RANGE/R_IN_0*8*(R_IN_0/8 + R_SOURCE)/128, //gain7 = x128
};
// ADC input range scale due to voltage dividers or preamplification on the board
#define ADC_0_CH0_SCALE_JO   0.004      //ADC#0 channel 0 with JP at I0 open or relay off, in volts
#define ADC_0_CH1_SCALE    (16000./39.*0.001)  //ADC#0 channel 1: High-voltage divider
#define ADC_1_CH0_SCALE    (10/2500.0)  //ADC#1 channel 0: uAmps
#define ADC_1_CH1_SCALE_JO   0.004      //ADC#1 channel 1 with JP5 open (voltage divider), in volts

// Definitions for AD5683 digital-to-analog converter !!!!!!!!!!!!!!!!  this stuff will go into the DAC library
#define AD5683_SPIMODE   1    //data accepted at falling edge of SCLK
#define AD5683_SET_DAC   0x30 //bits of highest (first) byte for writing & setting DAC
                              //lower 4 bits must be highest 4 bits of data
#define AD5683_RESET_MSB 0x40 //bits of highest (first) byte for reset (internal reference)

// Which hardware and closed jumpers do we have? Bits of present/closed jumpers are 1
#define ADC_0_PRESENT    0x01 //bit is set if ADC #0 was detected
#define ADC_1_PRESENT    0x02 //bit is set if ADC #1 was detected
#define LM35_PRESENT     0x04 //bit is set if LM35 temperature sensor was detected
#define RELAY_PRESENT    0x08 //bit is set if the relay for I0 input 2.5V/10V is present
#define JP_I0_CLOSED     0x10 //bit is set if JP3 7-8 is closed or relay on (to indicate 2.5V I0 range)
#define JP_AUX_CLOSED    0x20 //bit is set if JP5 is closed (to indicate 2.5V AUX range)

//SPI communication settings
SPISettings AD5683_SPI_SETTING(2500000, MSBFIRST, AD5683_SPIMODE); 
SPISettings AD7705_SPI_SETTING(2500000, MSBFIRST, AD7705_SPIMODE);

/* ---------- GLOBAL VARIABLES ---------- */ // !!!!!!!!!!!!!!! Will stay in main part, but see comments below
byte     hardwareDetected = 0;         //bits set indicate this hardware is spresent/jumper closed
byte     adc0Channel = AD7705_CH0;    //channel of ADC#0
byte     adc1Channel = AD7705_CH0;
byte     adc0Gain = 0;                //gain of ADC#0, 0...7
byte     adc1Gain = 0;
byte     adcUpdateRate = AD7705_50HZ; //update rate for both ADCs (will be set for line frequency) !!!!!!!!!!!!!!!!!! This is never used again, but probably should
uint16_t adc0RipplePP = 0;            //ripple (peak-peak) measured for ADC#0 during autogain at gain=0 !!!!!!!!!!!!!!!!!! Nor is this or the following. Not sure what Michael wanted to use these for. Perhaps as a substitute of summedMeasurements used only during the autogain?
uint16_t adc1RipplePP = 0;
bool adc0ShouldDecreaseGain = false;  //whether ADC#0 should increase its gain in the next cycle !!!!!!!!!!!!!!!!! Unused, same for the next one, but we should use them. However, right now we're just DECREASING the gain, not increasing it [except in findOptimalADCGains()]. Discuss with Michael.
bool adc1ShouldDecreaseGain = false;
// Self-calibration results are stored here
int32_t  selfCalDataForMedian[3][2][2]; //three values for median, two ADCs, last index is offset(0) & gain(1)
int32_t  selfCalDataVsGain[AD7705_MAX_GAIN+1][2][2][2]; //for each gain, two ADCs, two channels each, and last index is offset(0)&gain(1)
// Measurements of ADC#0, ADC#1 and LM35 sensor temperature are summed up here
int32_t summedMeasurements[N_MAX_MEAS];
int16_t measurement[N_MAX_MEAS];
//Union for Variables to convert from float to single bytes
union floatOrBytes{                  
  float asFloat; 
  byte asBytes[4];
} fDataOutput[N_MAX_MEAS];
//Variables to determing peak-peak distance in autogain
int16_t maximumPeak[2];
int16_t minimumPeak[2];
/* ---------- INITIALIZATION ---------- */
/** The setup routine runs once on power-up or on hardware reset */
void setup() {
  #ifdef DEBUG
    delay(1000);  // initialize serial communication (emulation on USB) at 9600 bits per second for debug
    Serial.begin(9600);
    delay(1000);
  #endif

  analogReference(INTERNAL);  // 2.56 V on Arduino Micro, for LM35 (and checking relay presence)
  //Set all inverted chip-select pins to high via pullup (to avoid low glitch)
  setChipSelectHigh(CS_DAC);
  setChipSelectHigh(CS_ADC_0);
  setChipSelectHigh(CS_ADC_1);
  pinMode(SCK, OUTPUT);  //should not be needed? but it did not work without
  pinMode(MOSI, OUTPUT);
  AD5683reset(CS_DAC);
  //Reset AD7705 IO. Note that (inverted) chip select = high does not reset it
  AD7705resetCommunication(CS_ADC_0);
  AD7705resetCommunication(CS_ADC_1);
  AD5683setVoltage(CS_DAC, 0x0000);
  #ifdef DEBUG
    Serial.print("hardware=0x");
    Serial.println(hardwareDetected, HEX);
  #endif
  #ifdef DEBUG
    pinMode(LED_BUILTIN, OUTPUT);
  #endif
}

void loop() {
/*
* Main loop with "readFromSerial()", which receives messages from 
* the computer and decodes it. "updateState()" interprets message
* and defines a state with defining the variable "currentState". Depending 
* on the state, the arduino continues with differnt functions
*/
  readFromSerial();
  updateState();

  switch (currentState){
    case STATE_SETUP_ADC:
      setUpAllADCs();
      break;
    case STATE_SET_VOLTAGE: 
      setVoltage(); 
      break;
    case STATE_TRIGGER_ADCS:
      if((millis()-currentTime) >= dac_settletime){
      triggerMeasurements();
      }
      break;
    case STATE_AUTOGAIN: 
      findOptimalADCGains(); 
      break; 
    case STATE_ADC_MEASURE:
      measureADCs();                                    // !!!!!!!!!!!!!! will probably be replaced in full by makeAndSumMeasurements() after edits (see comments further below)
      break;
    case STATE_ADC_VALUE_READY:
      sendAdcValue();
      break;
    case STATE_GET_HARDWARE:
       hardwareDetected = getHardwarePresent();
       encodeAndSend(hardwareDetected);
       initialCalibration();
      break;
    case STATE_ERROR:
      // !!!!!!!!! need to set up communication for this
      break;
   }
}

//================

void readFromSerial() {
/*
 * Receive one byte (i.e., character) from the serial line (i.e., PC) and store
 * it into "temp_buffer[]". Also track whether all bytes that need to be
 * received have been received by looking at whether the last byte read is
 * ENDMARKER. When this happens, temp_buffer will be:              !!!!!!!!!!!!!!!! Probably it would also be good to check here if the number of bytes read is the one expected?
 *     [STARTMARKER, byte with length of message, message, ENDMARKER]
 * 
 * Return into globals: 
 * -------------------
 * data_received_n: integer with length of message
 * temp_buffer: byte array
 */
  if(Serial.available() > 0) {                  // !!!!!!!!!!!!!!< this could become if(Serial.available()) [clearer]
    if(Serial.available() >= 64){Serial.println("Serial overflow."); currentState = STATE_ERROR;}
    byte byteRead = Serial.read();
    if (byteRead == STARTMARKER) {
      numBytesRead = 0; 
      readingFromSerial = true;
    }
    if(readingFromSerial) {
      inputRegister[numBytesRead] = byteRead;
      numBytesRead ++;
    } 
    if (byteRead == ENDMARKER) {
      readingFromSerial = false;
      newMessage = true;  // !!!!!! make sure to not always set to true
      msgLength = inputRegister[1]; 
      decodeMessage();
    }    
  }
}

void decodeMessage(){
/*
 * Puts received message from array "temp_buffer[]" into array 
 * "data_received[]", keeping only the actual characters (i.e., skipping
 * STARTMARKER == temp_buffer[0], the no.of bytes == temp_buffer[1], and 
 * ENDMARKER == temp_buffer[last]). The only special decoding is for bytes
 * that are SPECIAL_BYTE. In this case, the actual character is
 * SPECIAL_BYTE + the next character.
 * 
 * Return into globals: 
 * --------------------
 * data_received[]: byte array
 * data_received_counter: integer
 */
  data_received_counter = 0; 
  for (byte nthByte = 2; nthByte < numBytesRead; nthByte++) {
    byte x = inputRegister[nthByte];
    if (x == SPECIAL_BYTE) {
       nthByte++;
       x += inputRegister[nthByte];
    }
    data_received[data_received_counter] = x;
    data_received_counter++;
  }
}

//============================
void encodeAndSend(byte singleByte){
/*
 * Prepares message before sending it to the PC. Puts single 
 * byte into an array and forwards the array to the "real" 
 * encode message
 * This overloaded function essentially prepares a single-byte-long message
 * to be actually encoded in the next function. Having the two with the same
 * names prevents the rest of the code from having to figure out which function
 * to call depending on whether the message is a single byte or a byte array
 * 
 * Parameters: 
 * -----------
 * singleByte : single byte
 * 
 * Return into function: 
 * -----------
 * NewVar[]: byte array with one value
 */
  byte byteArray[1];
  byteArray[0] = singleByte;
  encodeAndSend(byteArray, 1);
}
void encodeAndSend(byte *byteArray, int len){
/*
 * Prepares message before sending it to the PC. Changes every 
 * byte which happens to have the same value as a STARTMARKER, an ENDMARKER or
 * a SPECIAL_BYTE to two bytes with a leading "SPECIAL_BYTE" and a following
 * "byte - SPECIAL_BYTE."
 * The message is put into the array "data_send[]". The length
 * of the array is defined in the variable "numBytesBeforeEncoding".
 * 
 * Parameters: 
 * -----------
 * byteArray : byte array
 * len : integer, length of byte array
 * 
 * Returns into globals
 * -----------
 * data_send: byte array
 */
  byte numBytesAfterEncoding = 0; 
  byte numBytesBeforeEncoding = 0; // !!!!!!! Can we remove this? Byte cast not clear
  for(int i = 0; i < len; i ++){
    if(byteArray[i] >= SPECIAL_BYTE){
      data_send[numBytesAfterEncoding] = SPECIAL_BYTE; 
      numBytesAfterEncoding++;
      data_send[numBytesAfterEncoding] = byteArray[i] - SPECIAL_BYTE; 
    }
    else{
      data_send[numBytesAfterEncoding] = byteArray[i]; 
    }
    numBytesAfterEncoding++;
  } 
  numBytesBeforeEncoding = len;
/*
 * Sends byte array "data_send" (i.e., the actual message) to PC as: 
 *   STARTMARKER
 *   numbers of bytes in actual message
 *   actual message
 *   ENDMARKER
 */
  Serial.write(STARTMARKER);
  Serial.write(numBytesBeforeEncoding);
  Serial.write(data_send, numBytesAfterEncoding);
  Serial.write(ENDMARKER);
}
//============================

void updateState() {
/*
 * Use the message received from the python script to  
 * update the state, provided that the Arduino is currently idle
 * 
 * Parameters from globals: 
 * ----------------
 * data_received[0]: byte array
 */
  if (newMessage){ // !!!!!!! need to remember to not always set newMessage = true
    switch(data_received[0]){
      case PC_HARDWARE:
        initialTime = millis(); 
        currentState = STATE_GET_HARDWARE;
        break;
      case PC_INIT_ADC:
        initialTime = millis();
        currentState = STATE_SETUP_ADC;
        break;
      case PC_SET_VOLTAGE:
        initialTime = millis();
        currentState = STATE_SET_VOLTAGE; 
        break;
      case PC_AUTOGAIN:
        initialTime = millis();
        resetMeasurementData();
        adc0Gain = 0;
        adc1Gain = 0;
        adcUpdateRate = AD7705_500HZ;
        selfCalibrateAllADCs(AD7705_500HZ);
        numMeasurementsToDo = 20; // !!!!!!!!! May need to be changed later on
        currentState = STATE_AUTOGAIN;       // !!!!!!!!!!!!!!!! something wrong here! measurement_n is still zero if the PC did not yet ask for a measurement. If it did, the number of measurements requested last time are used for the auto-gain. After a bit of fiddling, I realized that this value is set with setUpAllADCs(). This means that one cannot run this (nor the next one) unless the ADC has been initialized. We need to do some checking, probably on the fact that measurement_n should be a positive number (this means that setUpAllADCs has run). If not, I'd send back an ERROR. reset() should set the number back to zero.
        break;
      case PC_MEASURE:                       // !!!!!!!!!!!!!!!!  bug? Also here numMeasurementsToDo is not updated from anywhere, nor read. 
        initialTime = millis(); 
        currentState = STATE_ADC_MEASURE;
        break;
      case PC_RESET:
        reset();
        break;                             // !!!!!!!!!!!!!!!!!!! add one more case for the "?" command (or whichever character we will use to indicate the request of hardware configuration. This will need to call the getHardwarePresent() and return appropriate info to the PC. !!!!!!!!!!!!!!!!! And add Error state
    } 
    newMessage = false;
  }
}

//=========================

void findOptimalADCGains(){
/*
 * Find the optimal gain for all ADCs.
 * (1) measure several values. Remain in AUTOGAIN_STATE until done.            !!!!!!!!!!!!!!!  BUG: numMeasurementsToDo not set!
 * (2) for each ADC, pick a gain G such that the average of the measurements
 * is at least 25% of the range of values that can be measured with G
 * 
 * From Bernhard                                                               !!!!!!!!!!!!!!!
 * Function to initialse the Autogain. First it makes several ADC
 * measurements cycles, adds the values up to a sum and then  
 * devides by the number of measurement cycles. Then it decides 
 * for the best gain to start with. It decides for a gain, where 
 * the measured value is greater than 25% of the measurement 
 * range at the first time. 
 */
  int16_t autogain_value0;
  int16_t autogain_value1;
  measureADCsAndPeaks(); // !!!!!!!!!!!!! not at 500 Hz !!!!!!!!!!!!!!!!
  if(numMeasurementsDone == numMeasurementsToDo){
    autogain_value0 = max(abs(maximumPeak[0]),abs(minimumPeak[0])) + (maximumPeak[0] - minimumPeak[0]);
    autogain_value1 = max(abs(maximumPeak[1]),abs(minimumPeak[1])) + (maximumPeak[1] - minimumPeak[1]);
    numMeasurementsDone = 0;// !!!!!!!!!!!!!!!!!!! we should rather zero this counter later, after the correct gain has been found, just using resetMeasurementData() 
    if(hardwareDetected & ADC_0_PRESENT){
      while(((autogain_value0 << (adc0Gain + 1)) < ADC_RANGE_THRESHOLD) && (adc0Gain < AD7705_MAX_GAIN)){
        adc0Gain++;
      }
    }
    if (hardwareDetected & ADC_1_PRESENT){
      while(((autogain_value1 << (adc1Gain + 1)) < ADC_RANGE_THRESHOLD) && (adc1Gain < AD7705_MAX_GAIN)){
        adc1Gain++;
      }
    }
    encodeAndSend(PC_OK);
    currentState = STATE_IDLE; 
  }
  if((millis() -  initialTime) > timeout){               // !!!!!!!!!!!!!!!! perhaps go to an ERROR state?
    debugToPC("Timeout, Autogain failed, Arduino turns back to IDLE state"); 
    currentState = STATE_IDLE;
  }
}

void measureADCsAndPeaks(){
  if (hardwareDetected & ADC_0_PRESENT)
        {measurement[0] = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
        if(measurement[0] > maximumPeak[0]){
          maximumPeak[0] = measurement[0];
          }
        if(measurement[0] < minimumPeak[0]){
          minimumPeak[0] = measurement[0];
          }
        }
    if (hardwareDetected & ADC_1_PRESENT)
        {measurement[1] = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        if(measurement[1] > maximumPeak[1]){
          maximumPeak[1] = measurement[1];
          }
        if(measurement[1] < minimumPeak[1]){
          minimumPeak[1] = measurement[1];
          }
        }
numMeasurementsDone++;
}

//=========================
void setUpAllADCs(){ 
/* 
 *  Initialize the ADC from the parameters stored in the first 7 bytes of
 *  "data_received[]". The bytes have the following meaning:          // << TODO: list what each byte is used for
 */
  if(newMessage){
    adc0Channel = data_received[0];
    adc1Channel = data_received[1];
    adcUpdateRate = data_received[2]; 
    numMeasurementsToDo = data_received[3] << 8 | data_received[4];
    //maximum_gain = data_received[5]; //!!!!!!!!!! remove from Phython
    //selfCalibrateAllADCs(adcUpdateRate);
    setAllADCgainsAndCalibration();
    delay(1);
    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;//Back to command mode
    newMessage = false;
  }
  if((millis() -  initialTime) > timeout){           // !!!!!!!!!! Perhaps go to ERROR state here?
    debugToPC("Timeout, no ini Values for ADC received!"); 
    currentState = STATE_IDLE; //Back to cmd mode
    newMessage = false; 
  }
}
//=========================
void setVoltage(){
/*
 * Ask the DAC to provide a new voltage, if new settings are available.
 *
 * The voltage to be set is taken from the first two bytes of
 * the "data_received[]" buffer, that make up a uint16_t
 * The next two bytes of the buffer are packed into 
 * the unsigned long integer "dac_settletime", which defines 
 * the time interval that the Arduino will wait before deeming
 * the output voltage stable.
 *
 * The Arduino will go into the filter synchronization state 
 * STATE_TRIGGER_ADCS after the new voltage is set, triggering 
 * measurements on the next loop iteration
 */
  if(newMessage){  
    uint16_t dac_value = data_received[0] << 8 | data_received[1];
    dac_settletime = data_received[2] << 8 | data_received[3];  // !!!!!!!!!!!!! I don't understand why we're using an unsigned LONG type, and then accept only two bytes. Two bytes seem enough, I don't see the reason for waiting longer than 65.535 seconds, so perhaps we can change this into a uint16_t?
    AD5683setVoltage(CS_DAC, dac_value);
    currentState = STATE_TRIGGER_ADCS;
    currentTime = millis();
    newMessage = false;
    resetMeasurementData(); 
    if(adc0ShouldDecreaseGain){            // !!!!!!!!!!!!!!!! This section needs to be adapted to check whether either of the ADCs needs its gain to be reduced at this iteration, in which case the relevant self-calibration should be called
      // During the last measurement, the value read by the ADC has reached the
      // upper part of the range. The gain needs to be decreased, and the ADC
      // requires recalibration
      adc0Gain--;
      setAllADCgainsAndCalibration();
      //AD7705selfCalibrate(CS_ADC_0, adc0Channel, adc0Gain, adcUpdateRate); // !!!!!!!!!!!!!!!!!!!!!! ALL OR ONLY ONE? MAYBE CHECK BOTH AND IF ONE NEEDS TO BE CALIBRATED WE DO BOTH
      //AD7705waitForCalibration(CS_ADC_0, adc0Channel);
      delay(1);
      adc0ShouldDecreaseGain = false; 
    }
    if(adc1ShouldDecreaseGain){
      adc1Gain--;
      setAllADCgainsAndCalibration();
      //AD7705selfCalibrate(CS_ADC_1, adc1Channel, adc1Gain, adcUpdateRate); // !!!!!!!!!!!!!!!!!!!!!! ALL OR ONLY ONE? MAYBE CHECK BOTH AND IF ONE NEEDS TO BE CALIBRATED WE DO BOTH
      //AD7705waitForCalibration(CS_ADC_1, adc1Channel);
      delay(1);
      adc0ShouldDecreaseGain = false; 
    }
  }
  if((millis() -  initialTime) > timeout){                   // !!!!!!!!!!!!! Perhaps we should have an additional ERROR state that can do the reporting to the PC, then back to IDLE?
    debugToPC("Timeout, no ini Values for DAC received!"); 
    currentState = STATE_IDLE;
    newMessage = false; 
  }
}

void measureADCs(){
/*
 * This function acquires measurements from the ADCs.
 * (1) Take a defined number of measurements (as per numMeasurementsToDo), that are
 *     added together in 'summedMeasurements[]'.
 * (2) When the number of measurement points reaches the requested amount, the
 *     state of the Arduino is changed to ADC_VALUE_READY_STATE, that will
 *     handle the communication of the measurement to the PC.
 */
  makeAndSumMeasurements(); //                                           !!!!!!!!!!!!!!! summed Measurements not yet divided !!!!!!!!!!!!!!
  if(numMeasurementsDone == numMeasurementsToDo){
   // adc_value.asFloat /= numMeasurementsToDo;            // !!!!!!!!!!!!!! I think we can leave this part to the STATE_ADC_VALUE_READY handler !!!!!!!!!!!!!! USE MICHEALS ARRAY !!!!!!!!!!!!
    currentState = STATE_ADC_VALUE_READY;
  }
  if((millis() -  initialTime) > 10000){                   // !!!!!!!!!!!!!!!! probably here go to ERROR state too. Also, it's probably good to use the "timeout" variable instead of 10000
    debugToPC("Timeout, ADC measurement Timeout!"); 
    currentState = STATE_IDLE;
    newMessage = false; 
  }
}
//=========================
void sendAdcValue(){  // !!!!!!!!!!! this is the only place where we need to know what we want to measure exactly [except for which channel of the AD7705s, which is already used in initialiseADC()]. We should probably make it such that it can return several values (from the different ADCs/LM35) one after the other, if the PC asked for them. We can rename it to reflect this, e.g., "measurementsToPC()" or similar. The safest way would be to have this function look at WORKINPROGRESS (or however we will call it), and send the stuff needed from the current data read. We should consider if we should use several loop() iterations (at most N_MAX_MEAS, one per each ADC) to report the data back (keeping track of what needs still to be sent back by editing WORKINPROGRESS each time we send something) or if it's safer to do it in one shot, effectively stalling reading from the serial line in the meantime.
/*
 * This function sends the value to the PC by packing the 4 byte float 
 * number into an byte array.
 */
  getFloatMeasurements();
  encodeAndSend(fDataOutput[0].asBytes, 4); // !!!!!!!!!!!! need to implement way to only send demanded data
  encodeAndSend(fDataOutput[1].asBytes, 4);
  encodeAndSend(fDataOutput[2].asBytes, 4);
  resetMeasurementData();
  currentState = STATE_IDLE;
}
//=========================
void reset(){
/*
 * Resets the Arduino and sets it to IDLE-STATE
 */
  currentState = STATE_IDLE;  
  resetMeasurementData();
  dac_settletime = 0;
  currentTime = 0;  
  adc0Gain = 0;
  adc1Gain = 0;
  numMeasurementsToDo = 0;
  adc0ShouldDecreaseGain = false; 
  adc1ShouldDecreaseGain = false;
  AD7705resetCommunication(CS_ADC_0);
  AD7705resetCommunication(CS_ADC_1);
  AD5683reset(CS_DAC);
}
//=========================
void debugToPC(const char *debugmsg){  // !!!!!!!!!!!!!!! Since we're actually using this to send error messages back, I would rename it accordingly. Probably we should also define some error codes to accompany the message itself.
/*
 * Sends a debugmessage to the PC by masking the message 
 * with a STARTMARKER +  DEBUGBYTE + the debug message +
 * ENDMARKER
 * 
 * Parameters
 * ----------
 * debugmessage: const char array with variable field declaration
 */
  byte nb = PC_DEBUG;
  Serial.write(STARTMARKER);
  Serial.write(nb);
  Serial.println(debugmsg);
  Serial.write(ENDMARKER);
}
//=========================
void debugToPC(byte debugbyte){
  byte nb = 0;
  Serial.write(STARTMARKER);
  Serial.write(nb);
  Serial.write(debugbyte);
  Serial.write(ENDMARKER);
}

/** Resets the summed measurement data to 0 */
void resetMeasurementData() {                               // !!!!!!!!!! This needs to be called whenever we are done sending measured values back to the PC
    for (int i = 0; i < LENGTH(summedMeasurements); i++)
        summedMeasurements[i] = 0;
    numMeasurementsDone = 0;
}

/** Triggers the ADCs to start the measurements now (i.e., AD7705 STATE_TRIGGER_ADCS) */
void triggerMeasurements() {
    if (hardwareDetected & ADC_0_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    if (hardwareDetected & ADC_1_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);
    currentState = STATE_IDLE;
}

/** Waits for finishing the measurements, then adds the results for ADC#0, ADC#1
 *  and the temperature of the LM35 sensor (if present) to the global
 *  summedMeasurements array. Increments the global numMeasurements counter.
 *  as an array of three signed 16-bit integers.
 *  Results for non-existing components are 0.
 *  Does not trigger the ADCs. Note that ADCs have to be triggered after changing the channel. */
 // !!!!!!!!!!!!!!!!!!!!!!! This function needs editing, as it should check for both AD7705 whether the
 //    value read is OK with respect to output scale. In which case, it should
 //    either set a flag (for each ADC) that will trigger (in the next loop) a
 //    change of gain, or (if the value goes out of range) send an error message
 //    The logics should be similar to the one that we currently have in measureADCs()
 //    simple function to do that (e.g., checkMeasurementInADCRange(...) or similar)? needs: value, adc gain, pointer to new-gain flag
void makeAndSumMeasurements() {
    if (hardwareDetected & LM35_PRESENT)
        {measurement[2] = analogReadMedian(LM35_PIN);
        summedMeasurements[2] += measurement[2];}   //this one first while we probably have to wait for the others
    if (hardwareDetected & ADC_0_PRESENT)
        {measurement[0] = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
        checkMeasurementInADCRange(CS_ADC_0, adc0Channel, adc0Gain, &adc0ShouldDecreaseGain, measurement[0]);
        summedMeasurements[0] += measurement[0];
        }
    if (hardwareDetected & ADC_1_PRESENT)
        {measurement[1] = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        checkMeasurementInADCRange(CS_ADC_1, adc1Channel, adc1Gain, &adc1ShouldDecreaseGain, measurement[1]);
        summedMeasurements[1] += measurement[1];
        }
    numMeasurementsDone++;
}

void checkMeasurementInADCRange(byte chipSelectPin, byte channel, byte gain, bool* adcShouldDecreaseGain, int16_t adcValue){
/*
 * (1) If the value measured at the current call is larger than the
 *     ADC_RANGE_SATURATION_THRESHOLD threshold, the value is still acceptable, 
 *     but the gain will need to be decreased the next time that the DAC voltage 
 *     is increased. This is signaled by setting the adcShouldDecreaseGain flag, 
 *     which will be processed during the next call to initialiseDAC()
 * (2) If instead the value measured is truly at saturation (full scale), the
 *     measurement itself is incorrect, and needs to be repeated from scratch
 *     after decreasing the gain.
 * (3) If the value is at full scale but the gain cannot be decreased, an ERROR 
 *     state will be returned.
 */
   if(abs(adcValue)>ADC_RANGE_THRESHOLD && (gain > 0) && !*adcShouldDecreaseGain){
    // the measured value is above the "saturation" threshold, but not yet at
    // true saturation, which would make the measured value completely wrong.
    // Defer the decrease of gain to the next time the DAC voltage will be
    // increased
    *adcShouldDecreaseGain = true;     
    }
    if(((adcValue^=8000) == ADC_POSITIVE_SATURATION) || ((adcValue^=8000) == ADC_NEGATIVE_SATURATION)){
      if(gain>0){
        gain--;
        setAllADCgainsAndCalibration();
        //AD7705selfCalibrate(chipSelectPin, channel, gain, adcUpdateRate);
        //AD7705waitForCalibration(chipSelectPin, channel);
        delay(1);
        resetMeasurementData();
        *adcShouldDecreaseGain = false; 
      }
      if(gain==0){
        debugToPC("ADC Overflow, maximum gain reached, no serious measurements possible...");
        currentState = STATE_IDLE; 
      }
    }
}

/** Converts the global summedMeasurements to float in physical units:
 *  Volts, Amps, degC. */
void getFloatMeasurements() {                       // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! probably it is best to have this function rather write to the globals 3-element array of shared_memory type rather than passing an array and returning it. This way we can call this function once, only when reporting stuff back to the PC
    if (hardwareDetected & ADC_0_PRESENT) {
        fDataOutput[0].asFloat = summedMeasurements[0] * voltsPerBit[adc0Gain] / numMeasurementsToDo;
        if (adc0Channel == AD7705_CH0) {  //ADC#0 channel 0: I0 input (volts)
            if ((hardwareDetected & JP_I0_CLOSED) == 0)  //0-10V with jumper open
                fDataOutput[0].asFloat *= ADC_0_CH0_SCALE_JO;
        }
        else                            //ADC#0 channel 1: high voltage (volts)
            fDataOutput[0].asFloat *= ADC_0_CH1_SCALE;
    }
    if (hardwareDetected & ADC_1_PRESENT) {
        fDataOutput[1].asFloat = summedMeasurements[1] * voltsPerBit[adc1Gain] / numMeasurementsToDo;
        if (adc1Channel == AD7705_CH0) {  //ADC#1 channel 0: I0 at biased sample (amps)
            fDataOutput[1].asFloat *= ADC_1_CH0_SCALE;
        }
        else                            //ADC#1 channel 1: AUX
            if ((hardwareDetected & JP_AUX_CLOSED) == 0)  //0-10V with jumper open
                fDataOutput[1].asFloat *= ADC_1_CH1_SCALE_JO;
    }
    if (hardwareDetected & LM35_PRESENT) //LM35: degrees C
        fDataOutput[2].asFloat = summedMeasurements[2] * ARDUINO_ADU_DEG_C / numMeasurementsToDo;
}

/** Simultaneous self-calibration for both ADCs (if present) */
void selfCalibrateAllADCs(byte updateRate) {
    unsigned long startMillis = millis();
    if (hardwareDetected & ADC_0_PRESENT)
        AD7705selfCalibrate(CS_ADC_0, adc0Channel, adc0Gain, updateRate);
    if (hardwareDetected & ADC_1_PRESENT)
        AD7705selfCalibrate(CS_ADC_1, adc1Channel, adc1Gain, updateRate);
    if (hardwareDetected & ADC_0_PRESENT)
        AD7705waitForCalibration(CS_ADC_0, adc0Channel);
    if (hardwareDetected & ADC_1_PRESENT)
        AD7705waitForCalibration(CS_ADC_1, adc1Channel);
    unsigned long endMillis = millis();
//#ifdef DEBUG
//    Serial.print("t_selfCal=");
//    Serial.println(endMillis - startMillis);
//#endif
    delay(1); //AD7705 needs about 100 us to process the self-calibration result! 1 ms is on the safe side
}

/** Stores the result of the last self calibration of both ADCs if present */
void storeAllSelfCalibrationResults(int32_t targetArray[2][2]) {
  for (int iADC=0; iADC<2; iADC++) {
    bool adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;
    if (hardwareDetected & adcPresentMask) {
      byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;
      byte channel = iADC==0 ? adc0Channel : adc1Channel;
      int32_t offs = AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET);
      int32_t gain = AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN);
      targetArray[iADC][0] = offs;
      targetArray[iADC][1] = gain;
    }
  }
}

/** Sets the ADC gains according to the global adc0Gain, adc1Gain variables and restores
 *  the previously stored calibration values from the selfCalDataVsGain array.
 *  Also triggers the ADC measurements */
void setAllADCgainsAndCalibration() {
  for (int iADC=0; iADC<2; iADC++) {
    byte adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;
    if (hardwareDetected & adcPresentMask) {
      byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;
      byte channel = iADC==0 ? adc0Channel : adc1Channel;
      byte gain = iADC==0 ? adc0Gain : adc1Gain;
      int32_t cOffs = selfCalDataVsGain[gain][iADC][channel][0];
      int32_t cGain = selfCalDataVsGain[gain][iADC][channel][1];
      AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET, cOffs);
      AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN, cGain);
//      #ifdef DEBUG
//        Serial.print("gain=");
//        Serial.println(gain);
//        Serial.print(" cOffs=");
//        Serial.print(AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET));
//        Serial.print(" cGain=");
//        Serial.println(AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN));
//      #endif
      AD7705setGainAndTrigger(chipSelect, channel, gain);
    }
  }
}

/** Sets a digital output used as a chip select signal
 *  from the default high-impedance state to high (=unselected),
 *  without a glitch to the low state */
void setChipSelectHigh(byte ioPin) {
    pinMode(ioPin, INPUT_PULLUP);
    digitalWrite(ioPin, HIGH);
    pinMode(ioPin, OUTPUT);
}

/** Returns a bit mask of which hardware was detected
 *  checks for hardware present: ADCs, LM35, relay, jumpers */
uint16_t getHardwarePresent() {
    int result = 0;
    //Check for ADCs
    delay(10);    //make sure that input lines have settled (10 millisec)
    byte adc0comm = AD7705readCommRegister(CS_ADC_0, AD7705_CH0); //0xff if only pullup, no ADC
    if (adc0comm != 0xff)
        result |= ADC_0_PRESENT;
    byte adc1comm = AD7705readCommRegister(CS_ADC_1, AD7705_CH0);
    if (adc1comm != 0xff)
        result |= ADC_1_PRESENT;
    //Check for LM35 temperature sensor: the analog voltage should be within
    //the ADC range and settle to a similar value after connecting to a pullup resistor
    //(the LM35 cannot sink more than about 1 uA, so the pullup will drive it high)
    pinMode(LM35_PIN, INPUT);
    delay(10);    //make sure the voltage has settled (10 millisec)
    int sensorValue1 = analogReadMedian(LM35_PIN);
    pinMode(LM35_PIN, INPUT_PULLUP);  //measure the voltage with internal pullup
    delay(10);    //apply pullup for 10 millisec
    pinMode(LM35_PIN, INPUT);         //reset for usual measurements
    delay(10);    //make sure the voltage has settled (10 millisec)
    int sensorValue2 = analogReadMedian(LM35_PIN);
    if (sensorValue1 > 0 && sensorValue1 < LM35_MAX_ADU &&
        sensorValue2 < LM35_MAX_ADU && abs(sensorValue2 - sensorValue1) < 10)
        result |= LM35_PRESENT;
    //Check for relay present: if the relay is mounted, it should also have an
    //external pullup that results in about 0.12 V at the pin, about 48 ADUs
    int sensorValue = analogReadMedian(RELAY_PIN);
    if (sensorValue > RELAY_MIN_ADU && sensorValue < RELAY_MIN_ADU)
        result |= RELAY_PRESENT;
    else {
        //Check jumper at JP3 indicating 2.5 V I0 range set by user (if no relay)
        pinMode(JP_I0_PIN, INPUT_PULLUP);
        delay(1);
        if (digitalRead(JP_I0_PIN) == 0)
            result |= JP_I0_CLOSED;
        pinMode(JP_I0_PIN, INPUT);      //pullup off, reduces power consuption
    }
    //Check jumper JP5 indicating 2.5 V AUX range set by user
    pinMode(JP_AUX_PIN, INPUT_PULLUP);
    delay(1);
    if (digitalRead(JP_AUX_PIN) == 0)
        result |= JP_AUX_CLOSED;
    pinMode(JP_AUX_PIN, INPUT);
    return result;
}

/*initial calibration after getting hardware*/
void initialCalibration(){
 //initial self-calibration: done in parallel for both ADCs (if present) (old version of the calibration)
 //adc0Gain = 0;             adc1Gain = 0;
 //adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;
 //selfCalibrateAllADCs(adcUpdateRate);
 //adc0Channel = AD7705_CH1; adc1Channel = AD7705_CH1;
 //selfCalibrateAllADCs(adcUpdateRate);
 //adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;
 //triggerMeasurements();  //set to default channels
  
  //initial self-calibration: done in parallel for both ADCs (if present), for all gains (new version of the calibration)
  for(int channel=0; channel<2; channel++){
    if(channel==0) {adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;}
    if(channel==1) {adc0Channel = AD7705_CH1; adc1Channel = AD7705_CH1;}
    for (int gain=0; gain<=AD7705_MAX_GAIN; gain++) {
      adc0Gain = gain;
      adc1Gain = gain;
      for (int i=0; i<3; i++) {  //3 values for 3-point median
        selfCalibrateAllADCs(AD7705_50HZ);
        storeAllSelfCalibrationResults(selfCalDataForMedian[i]);
        #ifdef DEBUG
          Serial.print("gain=");
          Serial.print(gain);
          Serial.print(" cOffs=");
          Serial.print(selfCalDataForMedian[i][0][0]);
          Serial.print(" cGain=");
          Serial.println(selfCalDataForMedian[i][0][1]);
        #endif
      }
      for (int iADC=0; iADC<2; iADC++) {
        for (int offsNgain=0; offsNgain<2; offsNgain++) {
          int32_t median = getMedian32(
              selfCalDataForMedian[0][iADC][offsNgain],
              selfCalDataForMedian[1][iADC][offsNgain],
              selfCalDataForMedian[2][iADC][offsNgain]);
          
          selfCalDataVsGain[gain][iADC][channel][offsNgain] = median;
        }
      }
    }
  }
  adc0Gain = 0; adc1Gain = 0;
}
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
 /*void AD7705setClock(byte chipSelectPin, byte updateRate) {
   AD7705startIO(chipSelectPin);
   SPI.transfer(AD7705_REG_CLOCK);
   SPI.transfer(AD7705_CLK | updateRate);
   AD7705endIO(chipSelectPin);
 }*/

 /** Sets the gain (0...7) and triggers the measurement by touching STATE_TRIGGER_ADCS
  *  for the given channel of the AD7705 with the given chip select pin.
  *  Also sets bipolar unbuffered mode and normal operation,
  *  and reads the data register to make sure any old data go away */
void AD7705setGainAndTrigger(byte chipSelectPin, byte channel, byte gain) {
    AD7705startIO(chipSelectPin);
    //SPI.transfer16(0xffff); //reset should not be necessary
    //SPI.transfer16(0xffff);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_STATE_TRIGGER_ADCS | gain << 3);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(gain << 3);
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0xffff);
    AD7705endIO(chipSelectPin);
}

/** Sets gain (0...7) and update rate (AD7705_50HZ, AD7705_60HZ, or AD7705_500HZ) and
 *  performs self-calibration of the AD7705 with the given chip select pin.
 *  Self-calibration is done for the given channel (the calibration of the other channel
 *  remeins untouched).
 *  Thereafter, one has to wait until self-calibration is done, by calling
 *  AD7705waitForCalibration, or wait for the first data (AD7705waitAndReadData).
 *  Triggering (STATE_TRIGGER_ADCS) during the self-calibration phase is not allowed. */
void AD7705selfCalibrate(byte chipSelectPin, byte channel, byte gain, byte updateRate) {
    AD7705startIO(chipSelectPin);
    SPI.transfer(AD7705_REG_CLOCK | channel);
    SPI.transfer(AD7705_CLK | updateRate);     //set clock and update rate
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_SELFCAL | gain << 3);    //start self-calibration
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0xffff);                       //read data register to ensure DRDY is off.
    AD7705endIO(chipSelectPin);
}

/** Returns the value in an AD7705 offset or gain calibration register.
 *  'theRegister' must be AD7705_REG_OFFSET or AD7705_REG_GAIN.
 *  This function can be called after AD7705selfCalibrate to revrieve the calibration result */
int32_t AD7705getCalibrationRegister(byte chipSelectPin, byte channel, byte theRegister) {
  AD7705startIO(chipSelectPin);
  SPI.transfer(theRegister | AD7705_READ_REG | channel);
  int32_t result = (int32_t)SPI.transfer16(0xffff) << 8; //read high 16 bits
  result |= SPI.transfer(0xff);              //read low 8 bits
  AD7705endIO(chipSelectPin);
  return result;
}

/** Sets the value in an AD7705 offset or gain calibration register.
 *  'theRegister' must be AD7705_REG_OFFSET or AD7705_REG_GAIN.
 *  This function can be called to restore the result of a previous self-calibration */
int32_t AD7705setCalibrationRegister(byte chipSelectPin, byte channel, byte theRegister, int32_t value) {
  AD7705startIO(chipSelectPin);
  SPI.transfer(theRegister | channel);
  SPI.transfer16(value>>8);                  //write high 16 bits
  SPI.transfer(value&0xff);                  //write low 8 bits
  AD7705endIO(chipSelectPin);
}

/** Waits until the AD7705 with the given chip select pin has finished self-calibration */
void AD7705waitForCalibration(byte chipSelectPin, byte channel) {
    AD7705startIO(chipSelectPin);
    do {
        delayMicroseconds(10);
        SPI.transfer(AD7705_REG_SETUP | AD7705_READ_REG | channel);
    } while (SPI.transfer(0x0) & AD7705_SELFCAL);  //read setup register and check for self-calibration
    AD7705endIO(chipSelectPin);
}

/** Waits for data and then reads the given channel of the AD7705 with the given chip select pin.
 *  Returns a signed 16-bit int with the signed value (N.B. we use bipolar mode only) */
int16_t AD7705waitAndReadData(byte chipSelectPin, byte channel) {
    while (AD7705readCommRegister(chipSelectPin, channel) & AD7705_DRDY) {//DRDY bit is 0 when ready
        delayMicroseconds(100);
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

/* ---------- DAC AD5683 FUNCTIONS ---------- */

/** Resets the AD5683 DAC */
void AD5683reset(byte chipSelectPin) {
    SPI.beginTransaction(AD5683_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
    SPI.transfer(AD5683_RESET_MSB);
    SPI.transfer16(0);
    digitalWrite(chipSelectPin, HIGH);
}

/** Sets the AD5683 DAC output voltage to an unsigned 16-bit value */
void AD5683setVoltage(byte chipSelectPin, uint16_t dacValue) {
    uint16_t lowBytes = dacValue << 4;
    byte hiByte = (dacValue >> 12) | AD5683_SET_DAC;
    SPI.beginTransaction(AD5683_SPI_SETTING);
    digitalWrite(chipSelectPin, LOW);
    SPI.transfer(hiByte);
    SPI.transfer16(lowBytes);
    digitalWrite(chipSelectPin, HIGH);
}

/* ---------- UTILITY FUNCTIONS ---------- */

/** Reads the Arduiono ADC for a given pin (A0...) three times and returns the median */
uint16_t analogReadMedian(byte pin) {
  uint16_t a0 = analogRead(pin);
  uint16_t a1 = analogRead(pin);
  uint16_t a2 = analogRead(pin);
  return getMedian16(a0, a1, a2);
}

/** Gets the median of three numbers */
uint16_t getMedian16(uint16_t a0, uint16_t a1, uint16_t a2) {
  uint16_t maximum = biggest16(a0, a1, a2);
  if (maximum == a0) return bigger16(a1, a2);
  if (maximum == a1) return bigger16(a0, a2);
  else return bigger16(a0, a1);
}

uint16_t bigger16(uint16_t a, uint16_t b) {
  return (a > b) ? a : b;
}

uint16_t biggest16(uint16_t a, uint16_t b, uint16_t c) {
  return bigger16(a, bigger16(b,c));
}

  /** Gets the median of three numbers */
int32_t getMedian32(int32_t a0, int32_t a1, int32_t a2) {
  int32_t maximum = biggest32(a0, a1, a2);
  if (maximum == a0) return bigger32(a1, a2);
  if (maximum == a1) return bigger32(a0, a2);
  else return bigger32(a0, a1);
}

int32_t bigger32(int32_t a, int32_t b) {
  return (a > b) ? a : b;
}

int32_t biggest32(int32_t a, int32_t b, int32_t c) {
  return bigger32(a, bigger32(b,c));
}