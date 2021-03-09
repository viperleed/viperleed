/*
Main LeedControl File
---------------------
Author: Bernhard Mayr, Michael Schmid, Florian Dörr
Date: 03.03.2021
---------------------
*/
//Libraries
#include <Arduino.h>
#include <SPI.h>

//File with Arduino settings
//#include "viper-ino.h"														<< What did you think we should put in there? All the headers of the functions?
//#include "viper-ino_control.h"	    										<< How about this? How is it different from the previous one?

#define DEBUG           true  //debug mode, writes to serial line, for use in serial monitor

//Useful macros
#define LENGTH(array)  (sizeof(array) / sizeof((array)[0]))
#ifndef NAN						     										 // << This is never used
  #define NAN        (1.0/0)   //floating-point not-a-number
#endif

//global variables for communication with Python module
#define STARTMARKER 254
#define ENDMARKER 255
#define SPECIAL_BYTE 253
#define MAX_MESSAGE 16		    // << perhaps nicer MAX_MSG_LENGTH?

//definition of messages for communication with PC
#define PING 3
#define INIT_ADC_BYTE 4			// << Maybe we could rename this to PY_INIT_ADC? See also following for suggestions. We can replace PY (=python) with PC if it's clearer, or remove it altogether.
#define OK_BYTE 5				// << ?PY_OK?
#define ADC_MEASURE_BYTE 6		// << ?PY_MEASURE?
#define INIT_DAC_BYTE 7			// << ?PY_INIT_DAC or, better, PY_SET_VOLTAGE?
#define INIT_AUTOGAIN_BYTE 8	// << ?PY_AUTOGAIN?
#define RESET_BYTE 68			// << ?PY_RESET? Also, 68 is the ascii code for "D". Perhaps we should use something more explicative, like "R" (82)
#define DEBUG_BYTE 0			// << ?PY_DEBUG?
								// << We will need one more command (asking for hardware configuration). Perhaps we can use "?" (==63)

// << THINGS TO DO EVERYWHERE (to conform with c++ style):
// - use "lower-camel-case" for variables (i.e., thisIsANewVariable) instead of "lower-snake-case" (i.e., this_is_a_new_variable)
// - use "upper-camel-case" for constants (i.e., THIS_IS_A_CONSTANT). I think most are already correct
// - use "lower-camel-case" for functions

//Variables for the communication protocol
byte bytes_received = 0;            //counter for received bytes										  << probably better "numBytesRead"
byte data_received_n = 0;           //numbers of bytes received in message, 2nd byte of received message  << probably better "msgLength". This is set in getSerialData but never used outside of it. Perhaps we should use it to check that indeed we got the data that we expected.
byte data_received_counter = 0;     //total number of bytes received								      << This is also used only locally in decodeMessage, and it's probably better to keep it only local. There should be minimal to no performance hit by doing this. Also, renaming it to something clearer would help
byte data_received[MAX_MESSAGE];    //Real message														  << Judging from the code, and if I didn't misinterpret anything, I think we can have only one buffer for the input/output to the PC, so we can replace "data_received" and "data_send" with a single "message" array
byte data_send[MAX_MESSAGE];  
byte temp_buffer[MAX_MESSAGE];																		//    << Maybe rename as inputRegister?
byte data_send_count = 0;           //the number of 'real' bytes to be sent to the PC					  << See suggestion on dataToPC. This also needs to get a better name
byte data_total_send = 0;           //the number of bytes to send to PC taking account of encoded bytes   << also this one can become just a local, with a better name

//Finite state machine 
int FSM = 0; 																 // << This must become "currentState" or similar. FSM is unreadable. Perhaps all the following ones need to be changed to "STATE_*" (e.g., STATE_INIT_ADC), or we could even drop the "STATE" part, if it does not give confusion. Also, perhaps it would be better to have them #define(d) rather than as const int? (save a little memory, 2 bytes * 7 names)
const int IDLE_STATE = 0; 
const int INIT_ADC_STATE = 1;												 // << Maybe more explicative: STATE_SETUP_ADC?
const int INIT_DAC_STATE = 2;												 // << Maybe more explicative: STATE_SET_ENERGY or STATE_SET_VOLTAGE
const int FSYNC = 3; 														 // << This needs a better name, like STATE_TRIGGER_ADCS
const int ADC_MEASURE_STATE = 4; 
const int ADC_VALUE_READY_STATE = 5;
const int AUTOGAIN_STATE = 6; 

//Flags 
boolean in_progress = false;            //Flag True when Arduino receives new data  	  << This needs to become clearer, e.g. "readingFromSerial"
boolean new_data = false;               //Flag True when all new data has arrived   	  << Perhaps "newMessage", or "messageReady"
boolean autogain_measure = true;        //True when measurements are taken for auto-gain  << Perhaps "runningAutogain"?
boolean adc_newgain = false;            //True when a new gain should be set			  << This may need to become a 2-element array of bools, since we have two ADCs? Perhaps also rename to "needNewADCGain". Eventually it may be possible to move this into the ADC library. See also adc0ShouldIncreaseGain/adc1ShouldIncreaseGain from Michael's part

//ADC settings
byte gain = 0;																	// << this must become a two-byte array or two bytes, since we have two ADCs. Also needs a better name. See also the stuff that comes from Michael's globals. That's how we should do it. Will move into the ADC library then. Similar considerations for all the ADC variables below. See Michael's solution below
byte chn; 
byte polar_mode; 
byte update_rate;
byte maximum_gain = 7;
//byte minimum_gain = 0;

//ADC measurement settings
uint16_t measurement_n; 														// << better name: "numMeasurementsToDo"; needs to be unified with Michael's numMeasurements
uint16_t measurement_increment;													// << better name: "numMeasurementsDone"

//Fraction of the ADC scale that we consider saturation.
//Used to decide if the gain needs decreasing
float adc_saturation = 0.5;														// << this may also become #define(d), or at least const, since we do not change it anywhere. Perhaps also better to use a more explicative name, e.g. ADC_RANGE_SATURATION_THRESHOLD or similar?

//union to quickly convert values from float (4 bytes in memory) into its
//byte representation, by either interpreting the 4 bytes at once (float_value)
//or reading 4 bytes into a 4-byte array
union shared_memory{															// << this name is not great. Perhaps floatOrBytes?
  float float_value;  															// << maybe nicer to call this "asFloat"
  struct{
    byte byte_array[4]; 
  }byte_struct;      															// << and this "asBytes". Also, I don't understand why we need an extra struct here. I think a simple 4-byte array would do
}adc_value;			 															// << It is good to read from both ADCs at the same time. In this case we would need two of these adc_value registers. Actually, three, including the LM35. We can have an N_MAX_MEAS-long array of these to store the measured values, and pass it as an input to getFloatMeasurements (or, rather, make getFloatMeasurements access it from the globals)
//union shared_memory adc_value;												// << remove

//Timers (defined in milliseconds)
unsigned long starttime; 			  											// << perhaps "initialTime", or "initialMillis"
unsigned long timeout = 5000;         											// << this can become a #define, or at least a const  
unsigned long dac_settletime = 100;												// << not sure why we use a long here, as we accept only 2 bytes (uint16_t) in initialiseDAC()
unsigned long actualtime;			  											// << perhaps "currentTime" or "currentMillis"

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
#define WORKINPROGRESS     7   //Python script will decide which pin      		// << this needs to be NOT a #define because it needs updating. It should be a byte. Also, probably good to give it a more descriptive name, like, whatToMeasure. with 8 bits, we can have a bit representation of the devices. E.g.: ADC0 if WORKINPROGRESS == 0x01 == 0b00000001, ADC1 if WORKINPROGRESS == 0x02 == 0b00000010, LM35 if WORKINPROGRESS == 0x04 == 0b00000100. And combinations can be obtained by bitwise OR. In principle we could also encode it better, such that this single byte can be used also to determine which channel of which ADC to measure: e.g.: MEASURE_ADC0_CH0 = 0x01, MEASURE_ADC0_CH1 = 0x02, MEASURE_ADC1_CH0 = 0x04, MEASURE_ADC1_CH1 = 0x08, MEASURE_LM35 = 0x0f. Obviously then one needs to handle somehow the cases in which the python asks for measuring two channels of the same ADC, which cannot be done at the same time.

//Arduino internal reference voltage, ADC maximum, etc
#define VREF_INTERNAL    2.56  //arduino micro internal ADC reference is 2.56 V
#define ARDUINO_ADC_MAX  0x3ff //10-bit ADC
#define ARDUINO_ADU_VOLTS (VREF_INTERNAL/ARDUINO_ADC_MAX)      //volts per ADU
#define ARDUINO_ADU_DEG_C (100*VREF_INTERNAL/ARDUINO_ADC_MAX)  //LM35 degC per ADU
#define RELAY_MIN_ADU     24   //relay has about 0.12V with pullup, ADC signal is larger than this
#define RELAY_MAX_ADU     96   //ADC signal of relay with pullup is less than this
#define LM35_MAX_ADU     320   //LM35 signal should be less than 0.8 V (80 degC)

// Definitions for AD7705 analog-to-digital converter							<< All this stuff will go into a separate ADC library
#define AD7705_SPIMODE     3   //data accepted at rising edge of SCLK
//AD7705 communication register
#define AD7705_READ_REG  0x08  //set bit for reading a register
#define AD7705_CH0       0x00  //set bit for addressing channel 0
#define AD7705_CH1       0x01  //set bit for addressing channel 1
#define AD7705_REG_COMM  0x00  //bits to set for accessing communication register
#define AD7705_REG_SETUP 0x10  //bits to set for accessing setup register
#define AD7705_REG_CLOCK 0x20  //bits to set for accessing clock register
#define AD7705_REG_DATA  0x30  //bits to set for accessing data register (16 bit)
#define AD7705_DRDY      0x80  //bit mask for Data-ReaDY bit in communication register
//AD7705 setup register
#define AD7705_FSYNC     0x01  //bit to set for Filter SYNChronization in setup register
#define AD7705_BUFMODE   0x02  //bit to set for buffered mode (not used by us)
#define AD7705_UNIPOLAR  0x04  //bit to set for unipolar mode (not used by us)
#define AD7705_SELFCAL   0x40  //bit to set for self-calibration in setup register
//AD7705 clock register
#define AD7705_CLK       0x0C  //bits to set for 4.9152 MHz crystal: CLK=1, CLKDIV=1
#define AD7705_50HZ      0x04  //bits to set for 50 Hz update rate and suppression
#define AD7705_60HZ      0x05  //bits to set for 60 Hz update rate and suppression
#define AD7705_500HZ     0x07  //bits to set for 500 Hz update rate
// Correction for the finite source impedance and ADC input impedance			<< Since this stuff depends on our circuit, it will not be in the ADC library
#define R_SOURCE       1300.0  //source resistance of our circuit at ADC inputs, 1.3 kOhm
#define R_IN_0     (1.0/(7e-12*38400)) //input resistance of ADC7705 in gain0 = x1
#define REF_OVER_RANGE  (2.5/32768)    //volts per bit at gain0, bipolar, input 0ohm   << Ask Michael: shouldn't we have the actual reference voltage here?
const float voltsPerBit[] = {
        REF_OVER_RANGE * R_IN_0 / (R_IN_0 + R_SOURCE),         		 //gain0 = x1
        REF_OVER_RANGE * R_IN_0 / 2 / (R_IN_0 / 2 + R_SOURCE) / 2,   //gain1 = x2 has half R_in_0
        REF_OVER_RANGE * R_IN_0 / 4 / (R_IN_0 / 4 + R_SOURCE) / 4,   //gain2 = x4 has 1/4 R_in_0
        REF_OVER_RANGE * R_IN_0 / 8 / (R_IN_0 / 8 + R_SOURCE) / 8,   //gain3 = x8 and up: 1/8 R_in_0
        REF_OVER_RANGE * R_IN_0 / 8 / (R_IN_0 / 8 + R_SOURCE) / 16,  //gain4 = x16
        REF_OVER_RANGE * R_IN_0 / 8 / (R_IN_0 / 8 + R_SOURCE) / 32,  //gain5 = x32
        REF_OVER_RANGE * R_IN_0 / 8 / (R_IN_0 / 8 + R_SOURCE) / 64,  //gain6 = x64
        REF_OVER_RANGE * R_IN_0 / 8 / (R_IN_0 / 8 + R_SOURCE) / 128, //gain7 = x128
};
// Scaling factors for ADC input due to voltage dividers or preamplification on the board
#define ADC_0_CH0_SCALE_JO 4.0           //ADC#0 channel 0 with JP at I0 open or relay off
#define ADC_0_CH1_SCALE    (16000./39.)  //ADC#0 channel 1: High-voltage divider
#define ADC_1_CH0_SCALE    (1.e-5/2.5)   //ADC#1 channel 0: Amps
#define ADC_1_CH1_SCALE_JO 4.0           //ADC#1 channel 1 with JP5 open (voltage divider)

// Definitions for AD5683 digital-to-analog converter							<< this stuff will go into the DAC library
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

/* ---------- GLOBAL VARIABLES ---------- */		// << Will stay in main part, but see comments below
byte     hardwareDetected = 0;        //bits set indicate this hardware is present/jumper closed
byte     adc0Channel = AD7705_CH0;    //channel of ADC#0
byte     adc1Channel = AD7705_CH0;
byte     adc0Gain = 0;                //gain of ADC#0, 0...7
byte     adc1Gain = 0;
byte     adcUpdateRate = AD7705_50HZ; //update rate for both ADCs (will be set for line frequency)  	 << This is never used again, but probably should
uint16_t adc0RipplePP = 0;            //ripple (peak-peak) measured for ADC#0 during autogain at gain=0  << Nor is this or the following. Not sure what Michael wanted to use these for. Perhaps as a substitute of summedMeasurements used only during the autogain?
uint16_t adc1RipplePP = 0;
bool adc0ShouldIncreaseGain = false;  //whether ADC#0 should increase its gain in the next cycle		 << Unused, same for the next one, but we should use them. However, right now we're just DECREASING the gain, not increasing it [except in initialiseAutogain()]. Discuss with Michael.
bool adc1ShouldIncreaseGain = false;
// Measurements of ADC#0, ADC#1 and LM35 sensor temperature are summed up here
uint16_t numMeasurements;
int32_t  summedMeasurements[N_MAX_MEAS];
//================

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
  AD7705resetIO(CS_ADC_0);														// << maybe rename to "AD7705resetCommunication"
  AD7705resetIO(CS_ADC_1);
  //check for hardware present: ADCs, LM35, relay, jumpers
  hardwareDetected = getHardwarePresent();
  //initial self-calibration: done in parallel for both ADCs (if present)
  adc0Gain = 0;             adc1Gain = 0;
  adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;
  selfCalibrateAllADCs(AD7705_50HZ); 											// << Perhaps was meant to be adcUpdateRate? 
  adc0Channel = AD7705_CH1; adc1Channel = AD7705_CH1;
  selfCalibrateAllADCs(AD7705_50HZ);											// << Perhaps was meant to be adcUpdateRate?
  adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;  //set to default channels
  triggerMeasurements();
  AD5683setVoltage(CS_DAC, 0x0000);
  #ifdef DEBUG
    Serial.print("hw=0x");
    Serial.println(hardwareDetected, HEX);
    pinMode(LED_BUILTIN, OUTPUT);
  #endif
}

void loop() {
/*
* Main loop with "getSerialData()", which receives messages from 
* the computer and decodes them. "processData()" interprets message
* and defines the new state by changing the variable "FSM". Depending 
* on the state, the Arduino continues with different functions
*/
  if(FSM != FSYNC){
    getSerialData();  															// << maybe rename to readFromSerial() or readFromUSB()?
  }
  
  processData();

  switch (FSM){
    case INIT_ADC_STATE:
      initialiseADC();
      break;
    case INIT_DAC_STATE:
      initialiseDAC(); 															// << maybe rename to setEnergy() or setVoltage()? more explicative
      break;
    case FSYNC:
      if((millis()-actualtime) >= dac_settletime){
        DoFSYNC();																// << see my comment at DoFSYNC()
      }
      break;
    case AUTOGAIN_STATE: 
      initialiseAutogain(); 													// << perhaps findOptimalADCGains()?
      break; 
    case ADC_MEASURE_STATE:
      measureADC();																// << will probably be replaced in full by makeAndSumMeasurements() after edits (see comments further below)
      break;
    case ADC_VALUE_READY_STATE:
      sendAdcValue();
      break;
   }
}

//================

void getSerialData() {  													 // << maybe rename to readFromSerial() or readFromUSB()?
/*
 * Receive one byte (i.e., character) from the serial line (i.e., PC) and store
 * it into "temp_buffer[]". Also track whether all bytes that need to be
 * received have been received by looking at whether the last byte read is
 * ENDMARKER. When this happens, temp_buffer will be:							<< Probably it would also be good to check here if the number of bytes read is the one expected?
 *     [STARTMARKER, byte with length of message, message, ENDMARKER]
 * 
 * Return into globals: 
 * -------------------
 * data_received_n: integer with length of message
 * temp_buffer: byte array
 */
  if(Serial.available() > 0) {    											 // << this could become if(Serial.available()) [clearer]
    byte x = Serial.read();		  											 // << perhaps rename x to something better, like lastChar or byteRead
    if (x == STARTMARKER) {
      bytes_received = 0; 
      in_progress = true;
    }
    if(in_progress) {
      temp_buffer[bytes_received] = x;
      bytes_received++;
    } 
    if (x == ENDMARKER) {
      in_progress = false;
      new_data = true;  
      data_received_n = temp_buffer[1]; 
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
  for (byte n = 2; n < bytes_received; n++) {
    byte x = temp_buffer[n];													// << perhaps rename to nthByte or similar?
    if (x == SPECIAL_BYTE) {
       n++;
       x = x + temp_buffer[n];													// << x += temp_buffer[n]
    }
    //data_received[data_received_counter] = temp_buffer[n];					// << remove
    data_received[data_received_counter] = x;
    data_received_counter++;
  }
}

//============================
void encodeMessage(byte var){  												 // << this name (and the following one, that should be the same) is very counterintuitive to me. The functions both encode and send the message back. Perhaps rename to "encodeAndSend" or similar?
/*
 * Prepares message before sending it to the PC. Puts single 
 * byte into an array and forwards the array to the "real" 
 * encode message
 * This overloaded function essentially prepares a single-byte-long message
 * to be actually encoded in the next function. Having the two with the same
 * names prevents the rest of the code from having to figure out which function
 * to call depending on whether the message is a single byte or a byte array
 * 
 * 
 * Parameters: 
 * -----------
 * var : single byte															<< "var" is impossible to understand. Rename to singleByte, or onlyOneByte
 * 
 * Return into function: 
 * -----------
 * NewVar[]: byte array with one value											<< "NewVar" is impossible to understand. Rename to byteArray
 */
  byte NewVar[1];
  NewVar[0] = var;
  encodeMessage(NewVar, 1);
}
void encodeMessage(byte *var, int len){
/*
 * Prepares message before sending it to the PC. Changes every 
 * byte which happens to have the same value as a STARTMARKER, an ENDMARKER or
 * a SPECIAL_BYTE to two bytes with a leading "SPECIAL_BYTE" and a following
 * "byte - SPECIAL_BYTE."
 * The message is put into the array "data_send[]". The length
 * of the array is defined in the variable "data_total_send".
 * 
 * Parameters: 
 * -----------
 * var : byte array																<< "var" is impossible to understand. Rename to byteArray, or even "_message" (notice underscore, to avoid conflict with the "message" register), or "msg"
 * len : integer, length of byte array
 * 
 * Returns into globals
 * -----------
 * data_send: byte array
 * data_total_send: integer
 */
  data_send_count = 0;
  data_total_send = 0;
  for(int i = 0; i < len; i ++){
    if(var[i] >= SPECIAL_BYTE){
      data_send[data_send_count] = SPECIAL_BYTE; 
      data_send_count++;
      data_send[data_send_count] = var[i] - SPECIAL_BYTE; 
    }
    else{
      data_send[data_send_count] = var[i]; 
    }
    data_send_count++;
  } 
  data_total_send = len; 
  dataToPC(); 
}
void dataToPC() {															 // << this function is used only in the previous one, so we can incorporate it there. This way, the data_send_count becomes just local, and we can get rid of it from the globals, perhaps.
/*
 * Sends byte array "data_send" (i.e., the actual message) to PC as: 
 *   STARTMARKER
 *   numbers of bytes in actual message
 *   actual message
 *   ENDMARKER
 */
    Serial.write(STARTMARKER);
    Serial.write(data_total_send);
    Serial.write(data_send, data_send_count);
    Serial.write(ENDMARKER);
}
//============================

void processData() {													     // << perhaps rename it to updateState(), as this is what the function does?
/*
 * Use the message received from the python script to  
 * update the state, provided that the Arduino is currently idle
 * 
 * Parameters from globals: 
 * ----------------
 * data_received[0]: byte
 */
  if ((FSM == IDLE_STATE) && new_data){ 
    switch(data_received[0]){
      case PING:
        encodeMessage(PING);
        break;
      case INIT_ADC_BYTE:														// << perhaps this should be named STATE_SETUP_ADCS? See my comment on initialiseADC()
        starttime = millis();
        FSM = INIT_ADC_STATE;
        break;
      case INIT_DAC_BYTE:
        starttime = millis();
        FSM = INIT_DAC_STATE; 
        break;
      case INIT_AUTOGAIN_BYTE:
        starttime = millis(); 
        autogain_measure = true;
        measurement_increment = 0;												// << this and the next line can be replaced with a single call to resetMeasurementData()
        adc_value.float_value = 0;
        FSM = AUTOGAIN_STATE;													// << something wrong here! measurement_n is still zero if the PC did not yet ask for a measurement. If it did, the number of measurements requested last time are used for the auto-gain. After a bit of fiddling, I realized that this value is set with initialiseADC(). This means that one cannot run this (nor the next one) unless the ADC has been initialized. We need to do some checking, probably on the fact that measurement_n should be a positive number (this means that initialiseADC has run). If not, I'd send back an ERROR. ResetStates() should set the number back to zero.
        break;
      case ADC_MEASURE_BYTE:													// << bug? Also here measurement_n is not updated from anywhere, nor read. 
        starttime = millis(); 
        FSM = ADC_MEASURE_STATE;
        break;
      case RESET_BYTE:
        ResetStates();
        break;																	// << add one more case for the "?" command (or whichever character we will use to indicate the request of hardware configuration. This will need to call the getHardwarePresent() and return appropriate info to the PC.
    } 
    new_data = false;
  }
}

//=========================
void DoFSYNC(){																	// << This just needs to be replaced completely with Michael's function triggerMeasurements()
/*
 * Synchronization of the ADCs is used as a trigger for starting
 * the measurement. We do not read the result of the first measurement
 * though 
 
 * From Bernhard
 * Filter synchronization is proceed on the ADC. This is 
 * done whenever the DAC sets a new Value. In the end 
 * of the procedure, the data register of the ADC is
 * read-out and thrown away. 
 */
 // ad7705.FSYNC(chn, polar_mode, gain);                                                       ACHTUNG! DURCH WAS ERSETZEN?
  FSM = IDLE_STATE;
}
//=========================

void initialiseAutogain(){														// << This whole function needs to be adapted for processing both ADCs at the same time (checking if they are available). Probably also good to rename it to findOptimalADCGains(), as this does not actually set the gain.
/*
 * Find the optimal gain for all ADCs.
 * (1) measure several values. Remain in AUTOGAIN_STATE until done.				<< BUG: measurement_n not set!
 * (2) for each ADC, pick a gain G such that the average of the measurements
 *     is at least 25% of the range of values that can be measured with G
 
 * Bernhard:
 * First it makes several ADC
 * measurements cycles, adds the values up to a sum and then  
 * divides by the number of measurement cycles. Then it decides 
 * for the best gain to start with. It decides for a gain, where 
 * the measured value is greater than 25% of the measurement 
 * range at the first time. 
 */
  uint16_t adc_readvalue;
  float autogain_value; 
  if(autogain_measure){
    adc_readvalue = AD7705waitAndReadData(WORKINPROGRESS, chn);					// << replace this, and the next two lines with a single call to makeAndSumMeasurements()? This measures the LM35 too (if present) that is not needed, but it should not be a big deal (it should be faster than the other ADCs).
    adc_value.float_value += int(adc_readvalue^0x8000);
    measurement_increment++;
    if(measurement_increment == measurement_n){									// << BUG? measurement_n is initialized to zero at the beginning by the compiler, but not set to a finite number in processData() when the PC requires an auto-gain (the python is not sending how many measurements it should take)
      adc_value.float_value = abs(adc_value.float_value/measurement_n/32767);	// << not sure why there is a 32767 = 2^15-1 instead of 32768=2^15 that we expect for bipolar mode. check with Michael
      autogain_measure = false; 
      measurement_increment = 0;												// << we should rather zero this counter later, after the correct gain has been found, just using resetMeasurementData()
      gain = 0; 
    }
  }
  if(!autogain_measure){														// I'm not sure why this is done over several loop() iterations rather than having a while here (it's at most 7 iterations. The first part, i.e., the measurement itself, is not repeated). Then, "autogain_measure = true;", "encodeMessage(OK_BYTE);", and "FSM = IDLE_STATE;" can just go after the loops are over.
    autogain_value = adc_value.float_value*(1 << gain); 
    if(autogain_value > 0.25){
      autogain_measure = true;
      encodeMessage(OK_BYTE);
      encodeMessage(1<<gain); //sends actual gain, for debugging, delete		// << remove
      //ad7705.init(chn, polar_mode,  gain, update_rate);						// << remove
      FSM = IDLE_STATE; 
    }
    else{
      gain++;
      if(gain == maximum_gain){
        autogain_measure = true;
        encodeMessage(OK_BYTE);
        encodeMessage(1<<gain); //sends actual gain, for debugging, delete
        //ad7705.init(chn, polar_mode, gain, update_rate); 						// << remove
        FSM = IDLE_STATE; 
      }
    }
  }
  if((millis() -  starttime) > timeout){										// << perhaps go to an ERROR state?
    debugToPC("Timeout, Autogain failed, Arduino turns back to IDLE state"); 
    FSM = IDLE_STATE;
  }
}

//=========================
void initialiseADC(){  															// << Perhaps actually better to initialize all ADCs at once, and rename this function to setUpAllADCs(). The function is basically already there there from Michael, i.e., selfCalibrateAllADCs(updateRate).
/* 
 *  Initialize the ADC from the parameters stored in the first 7 bytes of
 *  "data_received[]". The bytes have the following meaning:					// << TODO: list what each byte is used for
 */
  if(new_data){
    chn = data_received[0];
    polar_mode = data_received[1]; 												// << Since we're not changing the polarity (it's already hard-coded into Michael's part to be bipolar) we can use this bit as the channel for the second ADC. Another good option would be to bit-encode the channel in the first byte, since we just have 2 ADCs. See my comment on WORKINPROGRESS above.
    update_rate = data_received[2]; 
    measurement_n = data_received[3] << 8 | data_received[4];
    maximum_gain = data_received[5];
    //encodeMessage(1<<gain);													// << remove
    AD7705selfCalibrate(WORKINPROGRESS, chn, gain, update_rate);				// << this needs to become a call to selfCalibrateAllADCs
    encodeMessage(OK_BYTE);
    FSM = IDLE_STATE;  //Back to command mode
    new_data = false;
  }
  if((millis() -  starttime) > timeout){
    debugToPC("Timeout, no ini Values for ADC received!"); 						// << Perhaps go to ERROR state here?
    FSM = IDLE_STATE; //Back to cmd mode
    new_data = false; 
  }
}
//=========================
void initialiseDAC(){															// << rename to setEnergy() or setVoltage()?
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
 * The Arduino will go into the filter synchronization state FSYNC 				<< rename!
 * after the new voltage is set, triggering measurements on the next
 * loop iteration
 */
  if(new_data){  
    uint16_t dac_value = data_received[0] << 8 | data_received[1];
    dac_settletime = data_received[2] << 8 | data_received[3];					// << I don't understand why we're using an unsigned LONG type, and then accept only two bytes. Two bytes seem enough, I don't see the reason for waiting longer than 65.535 seconds, so perhaps we can change this into a uint16_t?
    AD5683setVoltage(CS_DAC, dac_value);
    FSM = FSYNC;
    actualtime = millis();
    new_data = false;
    adc_value.float_value = 0;													// << replace this line and the next with a single call to resetMeasurementData()
    measurement_increment = 0; 
    if(adc_newgain){															// << This section needs to be adapted to check whether either of the ADCs needs its gain to be reduced at this iteration, in which case the relevant self-calibration should be called
	  // During the last measurement, the value read by the ADC has reached the
	  // upper part of the range. The gain needs to be decreased, and the ADC
	  // requires recalibration
      gain--;
      AD7705selfCalibrate(WORKINPROGRESS, chn, gain, update_rate);
      adc_newgain = false;
    }
  }
  if((millis() -  starttime) > timeout){ 										// << Perhaps we should have an additional ERROR state that can do the reporting to the PC, then back to IDLE?
    debugToPC("Timeout, no ini Values for DAC received!"); 
    FSM = IDLE_STATE;
    new_data = false; 
  }
}

void measureADC(){																// << This function needs to take measurements from all ADCs, maybe rename to "measureADCs"? see also my comments below
/*
 * This function acquires measurements from the ADC.
 * (1) Take a defined number of measurements (as per measurement_n), that are
 *     added together in 'adc_value.float_value'.
 * (2.1) If the value measured at the current call is larger than the
 *     adc_saturation threshold, the value is still acceptable, but the gain
 *     will need to be decreased the next time that the DAC voltage is
 *     increased. This is signaled by setting the adc_newgain flag, which will 
 *     be processed during the next call to initialiseDAC()
 * (2.2) If instead the value measured is truly at saturation (full scale), the
 *     measurement itself is incorrect, and needs to be repeated from scratch
 *     after decreasing the gain.
 * (3) When the number of measurement points reaches the requested amount, the
 *     state of the Arduino is changed to ADC_VALUE_READY_STATE, that will
 *     handle the communication of the measurement to the PC.
 */
/*
 * uint16_t adc_readvalue;
  float autogain_value; 
  if(autogain_measure){
    adc_readvalue = ad7705.readADResult(chn);
    Serial.println(adc_readvalue);//Debugging
    adc_value.float_value += int(adc_readvalue^0x8000);
    encodeMessage(adc_value.byte_struct.byte_array, 4);//Debugging
 */
  uint16_t adc_readvalue = AD7705waitAndReadData(WORKINPROGRESS, chn);			// << This line and the following three need to be replaced with a single call to makeAndSumMeasurements(). See also my comment on makeAndSumMeasurements(), as that function needs to include already the part that figures out whether we need a gain switch. This means that, actually, this whole function will be more or less replaced completely by makeAndSumMeasurements()
  int16_t adc_bipolarval = int(adc_readvalue^0x8000);							// << this is already done in AD7705waitAndReadData
  adc_value.float_value += adc_bipolarval;
  measurement_increment++;
  if(abs(adc_bipolarval)>(int16_t)(adc_saturation * 32767) && (gain > 0) && !adc_newgain){
	// the measured value is above the "saturation" threshold, but not yet at
	// true saturation, which would make the measured value completely wrong.
	// Defer the decrease of gain to the next time the DAC voltage will be
	// increased
    adc_newgain = true;         
  } 
  if(abs(adc_bipolarval) == (int16_t)32768 && gain == 0){									// << some notes: 32767 = 0x7FFF = 0xFFFF ^ 0x8000; 32768 = 0x8000
    debugToPC("ADC Overflow, minimum gain reached, no serious measurements possible...");	// << here we can go to an ERROR state
    FSM = IDLE_STATE; 
  }
  if(abs(adc_readvalue) == 0xFFFF && gain > 0){
	// the measured value is truly at saturation already. We cannot wait the
	// next voltage change to reduce the gain, as the measurement is already
	// wrong. The measurement needs to be repeated.
    gain--;
    AD7705selfCalibrate(WORKINPROGRESS, chn, gain, update_rate);
    adc_value.float_value = 0;													// << since we will have to repeat the measurements anyway, we can repeat them for all ADCs. This line and the next can be replaced with a single call to resetMeasurementData()
    measurement_increment = 0;
    adc_newgain = false; 
  }
  if(measurement_increment == measurement_n){
    adc_value.float_value /= measurement_n;										// << I think we can leave this part to the ADC_VALUE_READY_STATE handler
    FSM = ADC_VALUE_READY_STATE;
  }
  if((millis() -  starttime) > 10000){ 
    debugToPC("Timeout, ADC measurement Timeout!"); 							// << probably here go to ERROR state too. Also, it's probably good to use the "timeout" variable instead of 10000
    FSM = IDLE_STATE;
    new_data = false; 
  }
}
//=========================
void sendAdcValue(){															// << this is the only place where we need to know what we want to measure exactly [except for which channel of the AD7705s, which is already used in initialiseADC()]. We should probably make it such that it can return several values (from the different ADCs/LM35) one after the other, if the PC asked for them. We can rename it to reflect this, e.g., "measurementsToPC()" or similar. The safest way would be to have this function look at WORKINPROGRESS (or however we will call it), and send the stuff needed from the current data read. We should consider if we should use several loop() iterations (at most N_MAX_MEAS, one per each ADC) to report the data back (keeping track of what needs still to be sent back by editing WORKINPROGRESS each time we send something) or if it's safer to do it in one shot, effectively stalling reading from the serial line in the meantime.
/*
 * This function sends the value to the PC by packing the 4 byte float 
 * number into an byte array.
 */
  adc_value.float_value = (adc_value.float_value/32768)/(1<<gain);				// << this needs to become a call to getFloatMeasurements(). See also comments there
  encodeMessage(1<<gain);														// << not sure we need this. Probably not.
  encodeMessage(adc_value.byte_struct.byte_array, 4);							// << this will need to be replaced with the right number of calls, depending on the value of WORKINPROGRESS (either one after the other, or as separate loop() iterations)
  adc_value.float_value = 0; 													// this and the next line need to be replaced with a single call to resetMeasurementData() (when this call will happen depends on whether we do the reporting in a single function call or if we spread it over few loop() iterations)
  measurement_increment = 0; 
  FSM = IDLE_STATE;
}
//=========================
void ResetStates(){																// << maybe just reset()?
/*
 * Resets the Arduino and sets it to IDLE-STATE
 */
  FSM = IDLE_STATE;  
  adc_value.float_value = 0;													// call resetMeasurementData()
  measurement_increment = 0; 
  dac_settletime = 0;
  actualtime = 0;
  gain = 0;
																				// << I would also reset measurement_n here, probably. We should most likely consider resetting other quantities as well. e.g., all those in the "Flags" section
  adc_newgain = false; 
  AD7705resetIO(WORKINPROGRESS); //												// << This needs to be called for both ADCs
}
//=========================
void debugToPC(const char *debugmsg){											// << Since we're actually using this to send error messages back, I would rename it accordingly. Probably we should also define some error codes to accompany the message itself.
/*
 * Sends a debug message to the PC by masking the message 
 * with a STARTMARKER +  DEBUGBYTE + the debug message +
 * ENDMARKER
 * 
 * Parameters
 * ----------
 * debugmessage: const char array with variable field declaration
 */
  byte nb = DEBUG_BYTE;
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

// ------- Michael's functions -------


/** Resets the summed measurement data to 0 */
void resetMeasurementData() {												 // << This needs to be called whenever we are done sending measured values back to the PC
    for (int i = 0; i < LENGTH(summedMeasurements); i++)
        summedMeasurements[i] = 0;
    numMeasurements = 0;
}

/** Triggers the ADCs to start the measurements now (i.e., AD7705 FSYNC) */
void triggerMeasurements() {
    if (hardwareDetected & ADC_0_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    if (hardwareDetected & ADC_1_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);
}

/** Waits for finishing the measurements, then adds the results for ADC#0, ADC#1
 *  and the temperature of the LM35 sensor (if present) to the global
 *  summedMeasurements array. Increments the global numMeasurements counter.
 *  as an array of three signed 16-bit integers.
 *  Results for non-existing components are 0.
 *  Does not trigger the ADCs. Note that ADCs have to be triggered after changing the channel. */
 // << This function needs editing, as it should check for both AD7705 whether the
 //    value read is OK with respect to output scale. In which case, it should
 //    either set a flag (for each ADC) that will trigger (in the next loop) a
 //    change of gain, or (if the value goes out of range) send an error message
 //    The logics should be similar to the one that we currently have in measureADC()
 //    simple function to do that (e.g., checkMeasurementInADCRange(...) or similar)? needs: value, adc gain, pointer to new-gain flag
void makeAndSumMeasurements() {
    if (hardwareDetected & LM35_PRESENT)
        summedMeasurements[2] = analogReadMedian(LM35_PIN);   //this one first while we probably have to wait for the others	<< bug? needs to be +=
    if (hardwareDetected & ADC_0_PRESENT)
        summedMeasurements[0] = AD7705waitAndReadData(CS_ADC_0, adc0Channel);	// << bug? needs to be +=
    if (hardwareDetected & ADC_1_PRESENT)
        summedMeasurements[1] = AD7705waitAndReadData(CS_ADC_1, adc1Channel);   // << bug? needs to be +=
    numMeasurements++;
}

/** Converts the global summedMeasurements to float in physical units:
 *  Volts, Amps, degC. */
void getFloatMeasurements(float fOutput[]) {									// << probably it is best to have this function rather write to the globals 3-element array of shared_memory type rather than passing an array and returning it. This way we can call this function once, only when reporting stuff back to the PC
    if (hardwareDetected & ADC_0_PRESENT) {
        fOutput[0] = summedMeasurements[0] * voltsPerBit[adc0Gain] / numMeasurements;
        if (adc0Channel == AD7705_CH0) {  //ADC#0 channel 0: I0 input (volts)
            if ((hardwareDetected & JP_I0_CLOSED) == 0)  //0-10V with jumper open
                fOutput[0] *= ADC_0_CH0_SCALE_JO;
        }
        else                            //ADC#0 channel 1: high voltage (volts)
            fOutput[0] *= ADC_0_CH1_SCALE;
    }
    if (hardwareDetected & ADC_1_PRESENT) {
        fOutput[1] = summedMeasurements[1] * voltsPerBit[adc1Gain] / numMeasurements;
        if (adc1Channel == AD7705_CH0) {  //ADC#1 channel 0: I0 at biased sample (amps)
            fOutput[1] *= ADC_1_CH0_SCALE;
        }
        else                            //ADC#1 channel 1: AUX
            if ((hardwareDetected & JP_AUX_CLOSED) == 0)  //0-10V with jumper open
                fOutput[1] *= ADC_1_CH1_SCALE_JO;
    }
    if (hardwareDetected & LM35_PRESENT) //LM35: degrees C
        fOutput[2] = summedMeasurements[2] * ARDUINO_ADU_DEG_C / numMeasurements;
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
#ifdef DEBUG
    Serial.print("t_selfCal=");
    Serial.println(endMillis - startMillis);
#endif
}

/** Sets a digital output used as a chip select signal
 *  from the default high-impedance state to high (=unselected),
 *  without a glitch to the low state */
void setChipSelectHigh(byte ioPin) {
    pinMode(ioPin, INPUT_PULLUP);
    digitalWrite(ioPin, HIGH);
    pinMode(ioPin, OUTPUT);
}

/** Returns a bit mask of which hardware was detected */
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
void AD7705resetIO(byte chipSelectPin) {
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

 /** Sets the gain (0...7) and triggers the measurement by touching FSYNC
  *  for the given channel of the AD7705 with the given chip select pin.
  *  Also sets bipolar unbuffered mode and normal operation,
  *  and reads the data register to make sure any old data go away */
void AD7705setGainAndTrigger(byte chipSelectPin, byte channel, byte gain) {
    AD7705startIO(chipSelectPin);
    //SPI.transfer16(0xffff); //reset should not be necessary
    //SPI.transfer16(0xffff);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_FSYNC | gain << 3);
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(gain << 3);
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0x0);
    AD7705endIO(chipSelectPin);
}

/** Sets gain (0...7) and update rate (AD7705_50HZ, AD7705_60HZ, or AD7705_500HZ) and
 *  performs self-calibration of the AD7705 with the given chip select pin.
 *  Self-calibration is done for the given channel (the calibration of the other channel
 *  remains untouched).
 *  Thereafter, one has to wait until self-calibration is done, by calling
 *  AD7705waitForCalibration, or wait for the first data (AD7705waitAndReadData).
 *  Triggering (FSYNC) during the self-calibration phase is not allowed. */
void AD7705selfCalibrate(byte chipSelectPin, byte channel, byte gain, byte updateRate) {
    AD7705startIO(chipSelectPin);
    SPI.transfer(AD7705_REG_CLOCK | channel);
    SPI.transfer(AD7705_CLK | updateRate);     //set clock and update rate
    SPI.transfer(AD7705_REG_SETUP | channel);
    SPI.transfer(AD7705_SELFCAL | gain << 3);    //start self-calibration
    SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
    SPI.transfer16(0x0);                       //read data register to ensure DRDY is off.
    AD7705endIO(chipSelectPin);
}

/** Waits until the AD7705 with the given chip select pin has finished self-calibration */
void AD7705waitForCalibration(byte chipSelectPin, byte channel) {
    AD7705startIO(chipSelectPin);
    do {
        delayMicroseconds(100);
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
    int result = SPI.transfer16(0x0);
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
    byte result = SPI.transfer(0x0);
    AD7705endIO(chipSelectPin);
    return result;
}

/* ---------- DAC AD5683 FUNCTIONS ---------- */  // will move into the DAC library

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

/** Reads the Arduino ADC for a given pin (A0...) three times and returns the median */
uint16_t analogReadMedian(byte pin) {
    uint16_t a0 = analogRead(pin);
    uint16_t a1 = analogRead(pin);
    uint16_t a2 = analogRead(pin);
    return getMedian(a0, a1, a2);
}

/** Gets the median of three numbers */
uint16_t getMedian(uint16_t a0, uint16_t a1, uint16_t a2) {
    uint16_t max = biggest(a0, a1, a2);
    if (max == a0) return bigger(a1, a2);
    if (max == a1) return bigger(a0, a2);
    else return bigger(a0, a1);
}

uint16_t bigger(uint16_t a, uint16_t b) {
    return (a > b) ? a : b;
}

uint16_t biggest(uint16_t a, uint16_t b, uint16_t c) {
    return bigger(a, bigger(b, c));
}
