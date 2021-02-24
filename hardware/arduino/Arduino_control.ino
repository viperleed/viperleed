/*
  Test_ReadAnalogVoltage
  Test for ViPErLEED
*/

#include <Arduino.h>
#include <SPI.h>

/* ---------- CONSTANTS ---------- */

#define DEBUG           true  //debug mode, writes to serial line, for use in serial monitor

//Langage definitions
#define LENGTH(array)  (sizeof(array) / sizeof((array)[0]))
#ifndef NAN
  #define NAN        (1.0/0)   //floating-point not-a-number
#endif

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
//AD7705 communication register
#define AD7705_READ_REG  0x08  //set bit for reading a register
#define AD7705_CH0       0x00  //set bit for adressing channel 0
#define AD7705_CH1       0x01  //set bit for adressing channel 1
#define AD7705_REG_COMM  0x00  //bits to set for accessing communication register
#define AD7705_REG_SETUP 0x10  //bits to set for accessing setup register
#define AD7705_REG_CLOCK 0x20  //bits to set for accessing clock register
#define AD7705_REG_DATA  0x30  //bits to set for accessing data register (16 bit)
#define AD7705_DRDY      0x80  //bit mask for Data Ready bit in communication register
//AD7705 setup register
#define AD7705_FSYNC     0x01  //bit to set for FSYNC in setup register
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
#define R_IN_0     (1.0/(7e-12*38400)) //input resistance of ADC7705 in gain0 = x1
#define REF_OVER_RANGE  (2.5/32768)    //volts per bit at gain0, bipolar, input 0ohm
const float voltsPerBit[] = {
        REF_OVER_RANGE*R_IN_0/(R_IN_0 + R_SOURCE),         //gain0 = x1
        REF_OVER_RANGE*R_IN_0/2/(R_IN_0/2 + R_SOURCE)/2,   //gain1 = x2 has half R_in_0
        REF_OVER_RANGE*R_IN_0/4/(R_IN_0/4 + R_SOURCE)/4,   //gain2 = x4 has 1/4 R_in_0
        REF_OVER_RANGE*R_IN_0/8/(R_IN_0/8 + R_SOURCE)/8,   //gain3 = x8 and up: 1/8 R_in_0
        REF_OVER_RANGE*R_IN_0/8/(R_IN_0/8 + R_SOURCE)/16,  //gain4 = x16
        REF_OVER_RANGE*R_IN_0/8/(R_IN_0/8 + R_SOURCE)/32,  //gain5 = x32
        REF_OVER_RANGE*R_IN_0/8/(R_IN_0/8 + R_SOURCE)/64,  //gain6 = x64
        REF_OVER_RANGE*R_IN_0/8/(R_IN_0/8 + R_SOURCE)/128, //gain7 = x128
};
// ADC input range scale due to voltage dividers or preamplification on the board
#define ADC_0_CH0_SCALE_JO 4.0           //ADC#0 channel 0 with JP at I0 open or relay off
#define ADC_0_CH1_SCALE    (16000./39.)  //ADC#0 channel 1: High-voltage divider
#define ADC_1_CH0_SCALE    (1.e-5/2.5)   //ADC#1 channel 0: Amps
#define ADC_1_CH1_SCALE_JO 4.0           //ADC#1 channel 1 with JP5 open (voltage divider)

// Definitions for AD5683 digital-to-analog converter
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

/* ---------- GLOBAL VARIABLES ---------- */
byte     hardwareDetected = 0;         //bits set indicate this hardware is spresent/jumper closed
byte     adc0Channel = AD7705_CH0;    //channel of ADC#0
byte     adc1Channel = AD7705_CH0;
byte     adc0Gain = 0;                //gain of ADC#0, 0...7
byte     adc1Gain = 0;
byte     adcUpdateRate = AD7705_50HZ; //update rate for both ADCs (will be set for line frequency)
uint16_t adc0RipplePP = 0;            //ripple (peak-peak) measured for ADC#0 during autogain at gain=0
uint16_t adc1RipplePP = 0;
bool adc0ShouldIncreaseGain = false;  //whether ADC#0 should increase its gain in the next cycle
bool adc1ShouldIncreaseGain = false;
// Measurements of ADC#0, ADC#1 and LM35 sensor temperature are summed up here
uint16_t numMeasurements;
int32_t  summedMeasurements[N_MAX_MEAS];

//global variables for communication with Python-module
#define STARTMARKER 254
#define ENDMARKER 255
#define SPECIAL_BYTE 253
#define MAX_MESSAGE 16
#define PING 3
#define INIT_ADC_BYTE 4
#define OK_BYTE 5
#define ADC_MEASURE_BYTE 6
#define INIT_DAC_BYTE 7
#define INIT_AUTOGAIN_BYTE 8
#define RESET_BYTE 68
#define DEBUG_BYTE 0

//Variables for the communication protocoll
byte bytes_received = 0;            //counter for received bytes
byte data_received_n = 0;           //numbers of bytes received in message, 2nd byte of received message
byte data_received_counter = 0;     //total number of bytes received
byte data_received[MAX_MESSAGE];    //Real message
byte data_send[MAX_MESSAGE];  
byte temp_buffer[MAX_MESSAGE];
byte data_send_count = 0;           //the number of 'real' bytes to be sent to the PC
byte data_total_send = 0;           //the number of bytes to send to PC taking account of encoded bytes

//Flags 
boolean in_progress = false;            //Flag True when arduino receives new data
boolean new_data = false;               //Fag True when all new data has arrived

//Finite state machine 
int FSM = 0; 
const int IDLE_STATE = 0; 
const int INIT_ADC_STATE = 1;
const int INIT_DAC_STATE = 2;
const int FSYNC = 3; 
const int ADC_MEASURE_STATE = 4; 
const int ADC_VALUE_READY_STATE = 5;
const int AUTOGAIN_STATE = 6; 

/* ---------- INITIALIZATION ---------- */
/** The setup routine runs once on power-up or reset */
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
  AD7705resetIO(CS_ADC_0);
  AD7705resetIO(CS_ADC_1);
  //check for hardware present: ADCs, LM35, relay, jumpers
  hardwareDetected = getHardwarePresent();
  //initial self-calibration: done in parallel for both ADCs (if present)
  adc0Gain = 0;             adc1Gain = 0;
  adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;
  selfCalibrateAllADCs(AD7705_50HZ);
  adc0Channel = AD7705_CH1; adc1Channel = AD7705_CH1;
  selfCalibrateAllADCs(AD7705_50HZ);
  adc0Channel = AD7705_CH0; adc1Channel = AD7705_CH0;
  triggerMeasurements();  //set to default channels
  AD5683setVoltage(CS_DAC, 0x0000);
  #ifdef DEBUG
    Serial.print("hw=0x");
    Serial.println(hardwareDetected, HEX);
    pinMode(LED_BUILTIN, OUTPUT);
  #endif
}

/** Main loop by Michael*/

/** The main loop routine runs over and over again forever 
void loop() {
  digitalWrite(LED_BUILTIN, HIGH);   // turn the LED on (HIGH is the voltage level)
  unsigned long startMillis = millis(); 
  resetMeasurementData();
  //triggerMeasurements();
  makeAndSumMeasurements();
  unsigned long endMillis = millis();
  float fData[N_MAX_MEAS] = {NAN, NAN, NAN};
  getFloatMeasurements(fData);
  Serial.print("dt=");
  Serial.print(endMillis - startMillis);
  for (byte i=0; i<N_MAX_MEAS; i++) {
    Serial.print(" val[");
    Serial.print(i);
    Serial.print("]=");
    byte nDigits = -log10(fData[i])+4;
    if (nDigits < 0) nDigits = 0;
    Serial.print(fData[i], nDigits);
  }
  Serial.println();
  digitalWrite(LED_BUILTIN, LOW);    // turn the LED off by making the voltage LOW
  delay(1000);                       // wait for a second
}*/

/**Main loop by Bernhard*/

void loop() {
/*
* Main loop with "getSerialData()", which receives messages from 
* the computer and decodes it. "processData()" interprets message
* and defines a state with defining the variable "FSM". Depending 
* on the state, the arduino continues with differnt functions
*/
  if(FSM != FSYNC){
    getSerialData();
  }
  
  processData();

  switch (FSM){
    case INIT_ADC_STATE:
      initialiseADC();
      break;
    case INIT_DAC_STATE: 
      initialiseDAC(); 
      break;
    case FSYNC:
      if((millis()-actualtime) >= dac_settletime){
        DoFSYNC();
      }
      break;
    case AUTOGAIN_STATE: 
      initialiseAutogain(); 
      break; 
    case ADC_MEASURE_STATE:
      measureADC();
      break;
    case ADC_VALUE_READY_STATE:
      sendAdcValue();
      break;
   }
}

/** Resets the summed measurement data to 0 */
void resetMeasurementData() {
  for (int i=0; i<LENGTH(summedMeasurements); i++)
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
void makeAndSumMeasurements() {
  if (hardwareDetected & LM35_PRESENT)
    summedMeasurements[2] = analogReadMedian(LM35_PIN);   //this one first while we probably have to wait for the others
  if (hardwareDetected & ADC_0_PRESENT)
    summedMeasurements[0] = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
  if (hardwareDetected & ADC_1_PRESENT)
    summedMeasurements[1] = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
  numMeasurements++;
}

/** Converts the global summedMeasurements to float in physical units:
 *  Volts, Amps, degC. */
void getFloatMeasurements(float fOutput[]) {
  if (hardwareDetected & ADC_0_PRESENT) {
    fOutput[0] = summedMeasurements[0] * voltsPerBit[adc0Gain] / numMeasurements;
    if (adc0Channel == AD7705_CH0) {  //ADC#0 channel 0: I0 input (volts)
      if ((hardwareDetected & JP_I0_CLOSED) == 0)  //0-10V with jumper open
        fOutput[0] *= ADC_0_CH0_SCALE_JO;
    } else                            //ADC#0 channel 1: high voltage (volts)
      fOutput[0] *= ADC_0_CH1_SCALE;
  }
  if (hardwareDetected & ADC_1_PRESENT) {
    fOutput[1] = summedMeasurements[1] * voltsPerBit[adc1Gain] / numMeasurements;
    if (adc1Channel == AD7705_CH0) {  //ADC#1 channel 0: I0 at biased sample (amps)
      fOutput[1] *= ADC_1_CH0_SCALE;
    } else                            //ADC#1 channel 1: AUX
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
  SPI.transfer(AD7705_FSYNC | gain<<3);
  SPI.transfer(AD7705_REG_SETUP | channel);
  SPI.transfer(gain<<3);
  SPI.transfer(AD7705_REG_DATA | AD7705_READ_REG | channel);
  SPI.transfer16(0x0);
  AD7705endIO(chipSelectPin);
}

/** Sets gain (0...7) and update rate (AD7705_50HZ, AD7705_60HZ, or AD7705_500HZ) and
 *  performs self-calibration of the AD7705 with the given chip select pin.
 *  Self-calibration is done for the given channel (the calibration of the other channel
 *  remeins untouched).
 *  Thereafter, one has to wait until self-calibration is done, by calling
 *  AD7705waitForCalibration, or wait for the first data (AD7705waitAndReadData).
 *  Triggering (FSYNC) during the self-calibration phase is not allowed. */
void AD7705selfCalibrate(byte chipSelectPin, byte channel, byte gain, byte updateRate) {
  AD7705startIO(chipSelectPin);
  SPI.transfer(AD7705_REG_CLOCK | channel);
  SPI.transfer(AD7705_CLK | updateRate);     //set clock and update rate
  SPI.transfer(AD7705_REG_SETUP | channel); 
  SPI.transfer(AD7705_SELFCAL | gain<<3);    //start self-calibration
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
  while(AD7705readCommRegister(chipSelectPin, channel) & AD7705_DRDY) {//DRDY bit is 0 when ready
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
  return bigger(a, bigger(b,c));
}

/* -------- BERNHARD'S FUNCTIONS -------- */

void getSerialData() {
/*
 * Function receives bytes from computer and packs whole message 
 * into the array "temp_buffer[]". Array consists then out of 
 * STARTMARKER + byte with length of message + message + 
 * ENDMARKER
 * 
 * Return into globals: 
 * -------------------
 * data_received_n: integer with length of message
 * temp_buffer: byte array
 */
  if(Serial.available() > 0) {    
    byte x = Serial.read();
    if (x == STARTMARKER) {
      bytes_received = 0; 
      in_progress = true;
    }
    if(in_progress) {
      temp_buffer[bytes_received] = x;
      bytes_received ++;
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
 * "data_received[]", but just the messsage (without STARTMARKER, 
 * ENDMARKER). It decodes also all the SPECIAL_BYTES to bytes 
 * which are STARTMARKER,ENDMARKER or SPECIAL_BYTES
 * 
 * Return into globals: 
 * --------------------
 * data_received[]: byte array
 * data_received_counter: integer
 */
  data_received_counter = 0; 
  for (byte n = 2; n < bytes_received; n++) {
    byte x = temp_buffer[n];
    if (x == SPECIAL_BYTE) {
       n++;
       x = x + temp_buffer[n];
    }
    //data_received[data_received_counter] = temp_buffer[n];
    data_received[data_received_counter] = x;
    data_received_counter++;
  }
}

void encodeMessage(byte var){
/*
 * Prepares message before sending it to the PC. Puts single 
 * byte into an array and forwards the array to the "real" 
 * encode message
 * 
 * Parameters: 
 * -----------
 * var : single byte
 * 
 * Return into function: 
 * -----------
 * NewVar[]: byte array with one value
 */
  byte NewVar[1];
  NewVar[0] = var;
  encodeMessage(NewVar, 1);
}

void encodeMessage(byte *var, int len){
/*
 * Prepares Message before sending it to the PC. Changes every 
 * byte which is a STARTBYTE, ENDBYTE or SPECIAL_BYTE to two 
 * bytes with a leading SPECIAL_BYTE and a byte - SPECIAL_BYTE.  
 * The Message is put into the array "data_send[]". The length
 * of the array is defined in the variable "data_total_send".
 * 
 * Parameters: 
 * -----------
 * var : byte array
 * len : integer, length of bytearray
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

void dataToPC() {
/*
 * Sends bytearray "data_send" to PC by masking message, 
 * beginning with STARTMARKER + numbers of bytes in message +
 * message + ENDMARKER
 */
    Serial.write(STARTMARKER);
    Serial.write(data_total_send);
    Serial.write(data_send, data_send_count);
    Serial.write(ENDMARKER);
}

void processData() {
/*
 * When the arduino is in IDLE_STATE, it turns the 
 * received message from "data_received[]" into 
 * a new state of the arduino.
 * 
 * Parameters from globals: 
 * ----------------
 * data_received[0]: byte array
 */
  if ((FSM == IDLE_STATE) && new_data){ 
    switch(data_received[0]){
      case PING:
        encodeMessage(PING);
        break;
      case INIT_ADC_BYTE:
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
        measurement_increment = 0;
        adc_value.float_value = 0;
        FSM = AUTOGAIN_STATE;
        break;
      case ADC_MEASURE_BYTE:
        starttime = millis(); 
        FSM = ADC_MEASURE_STATE;
        break;
      case RESET_BYTE:
        ResetStates();
        break;
    } 
    new_data = false;
  }
}

void initialiseADC(){ 
/* 
 *  Function initialise the ADC. First it takes the received message 
 *  "data_received[]" and puts all the bytes into a property of 
 *  the ADC. The message must consist out of 7 bytes. Then it 
 *  sends the initialise command to the ADC
 */
  if(new_data){
    chn = data_received[0];
    polar_mode = data_received[1]; 
    update_rate = data_received[2]; 
    measurement_n = data_received[3] << 8 | data_received[4];
    maximum_gain = data_received[5];
    //encodeMessage(1<<gain);
   // ad7705.init(chn, polar_mode, gain, update_rate);                                       ERSETZEN!!!!!!!!!!!!!!!!!!!
   //AD7705selfCalibrate(byte chipSelectPin, chn, gain, update_rate)                 DURCH DAS/MUSS ADC AUSWÄHLEN KÖNNEN
    encodeMessage(OK_BYTE);
    FSM = IDLE_STATE;//Back to command mode
    new_data = false;
  }
  if((millis() -  starttime) > timeout){
    debugToPC("Timeout, no ini Values for ADC received!"); 
    FSM = IDLE_STATE; //Back to cmd mode
    new_data = false; 
  }
}

void ResetStates(){
/*
 * Resets the Arduino and sets it to IDLE-STATE
 */
  FSM = IDLE_STATE;  
  adc_value.float_value = 0;
  measurement_increment = 0; 
  dac_settletime = 0;
  actualtime = 0;  
  gain = 0;
  adc_newgain = false; 
  ad7705.reset();
}
//=========================
void debugToPC(const char *debugmsg){
/*
 * Sends a debugmessage to the PC by masking the message 
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
