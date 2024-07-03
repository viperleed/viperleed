/*
ViPErLEED - Driver for INA239
---------------------
Author: Stefan Mitterhöfer
Date: 26.04.2024
---------------------
*/

#ifndef _VIPERLEED_INA239
#define _VIPERLEED_INA239

#include <SPI.h>
#define DEBUG false

#define INA_1_SPI_CS        7   // Pin for SPI communication; PE6 on ATmega32U4
#define INA_2_SPI_CS        8   // Pin for SPI communication; PB4 on ATmega32U4

/*
Communication parameters for the INA239:
  - The INA239 always operates in slave mode (read-only)
  - Baud rate: 10 MBaud/s max.
  - MSbit first, clock polarity (CPOL=0) and phase (CPHA=1)
*/
#define INA239_SPI_BAUD  1E7            // Can be increased up to 10 MHz
#define INA239_SPIMODE   SPI_MODE1      // SPI_MODE1 means CPOL=0, CPHA=1
#define INA239_SPI_SETTINGS SPISettings(INA239_SPI_BAUD, MSBFIRST, INA239_SPIMODE)

// INA239 constants
#define INA239_I_MAX        5           // Maximum expected Current in Amperes
#define INA239_R_SHUNT      0.01        // Resistance in Ohm
// Shunt full scale range selection across IN+ (PIN10) and IN– (PIN9).
// 0 = ±163.84 mV, 1 = ±40.96 mV, see datasheet page 20.
#define INA239_ADC_RANGE    1

// Resolution of the VBUS register, Conversion factor: 3.125 mV/LSB,
// see datasheet page 22. 
#define INA239_Voltage_LSB      0.003125
// Resolution of the DIETEMP register, Conversion factor: 125 m°C/LSB,
// see datasheet page 23. 
#define INA239_Temperature_LSB  0.125
// Resolution of the CURRENT register, see datasheet pages 27f.
#define INA239_Current_LSB      (INA239_I_MAX/pow(2, 15))

// SHUNT_CAL provides the device with a conversion constant value that 
// represents shunt resistance used to calculate current value in Amperes
// see datasheet pages 27f.
#define INA239_SHUNT_CAL_0  (819.2*pow(10, 6)*INA239_Current_LSB*INA239_R_SHUNT)

#if INA239_ADC_RANGE == 0
  // Shunt voltage conversion factor: 5 µV/LSB
  #define INA239_VSHUNT_LSB   0.005
  #define INA239_SHUNT_CAL    INA239_SHUNT_CAL_0
#elif INA239_ADC_RANGE == 1
  // Shunt voltage conversion factor: 1.25 µV/LSB
  #define INA239_VSHUNT_LSB   0.00125
  // SHUNT_CAL must be multiplied by 4 for ADCRANGE = 1
  #define INA239_SHUNT_CAL    (INA239_SHUNT_CAL_0*4)
#endif

/* 
// INA239 SPI Register addresses, see datasheet page 19.
// For detailed desciption see datasheet pages 20ff.
CONFIG              0x00  // General configuration
ADC_CONFIG          0x01  // Specific configuration of the measurement parameters
SHUNT_CAL           0x02  // Conversion constant value of shunt resistance 
VSHUNT              0x04  // Differential voltage measured across the shunt in mV
VBUS                0x05  // Voltage at PIN8 in Volts
DIETEMP             0x06  // Internal device temperature measurement in °C
CURRENT             0x07  // Calculated current output in Amperes
POWER               0x08  // NOT USED; Calculated power output 
DIAG_ALERT          0x0B  // NOT USED; Various diagnostics and alerts, see pages 23ff
SHUNT_OVERVOLTAGE   0x0C  // NOT USED; Shunt Overvoltage Threshold
SHUNT_UNDERVOLTAGE  0x0D  // NOT USED; Shunt Undervoltage Threshold
BUS_OVERVOLTAGE     0x0E  // NOT USED; PIN8 Overvoltage Threshold
BUS_UNDERVOLTAGE    0x0F  // NOT USED; PIN8 Undervoltage Threshold
TEMP_LIMIT          0x10  // NOT USED; Temperature Over-Limit Threshold
POWER_LIMIT         0x11  // NOT USED; Power Over-Limit Threshold
MANUFACTURER_ID     0x3E  // NOT USED; Manufacturer ID
DEVICE_ID           0x3F  // NOT USED; Device identification and revision
*/

// INA239 SPI Frames (8-bit)  READ:  ADDR5 ADDR4 ADDR3 ADDR2 ADDR1 ADDR0 0 1
// see datasheet pages 18f.
#define INA239_READ_CONFIG          (0x00 << 2 | 0x01)
#define INA239_READ_ADC_CONFIG      (0x01 << 2 | 0x01)
#define INA239_READ_SHUNT_CAL       (0x02 << 2 | 0x01)
#define INA239_READ_VSHUNT          (0x04 << 2 | 0x01)
#define INA239_READ_VOLTAGE         (0x05 << 2 | 0x01)
#define INA239_READ_TEMPERATURE     (0x06 << 2 | 0x01)
#define INA239_READ_CURRENT         (0x07 << 2 | 0x01)

// INA239 SPI Frames (8-bit) WRITE: ADDR5 ADDR4 ADDR3 ADDR2 ADDR1 ADDR0 0 0
// see datasheet pages 18f.
#define INA239_WRITE_CONFIG         (0x00 << 2 | 0x00)
#define INA239_WRITE_ADC_CONFIG     (0x01 << 2 | 0x00)
#define INA239_WRITE_SHUNT_CAL      (0x02 << 2 | 0x00)

// INA239 SPI Commands (16-bit)
#define INA239_CONFIG_RESET           (0x01 << 15) // Set bit 15 HIGH in Register CONFIG
#if INA239_ADC_RANGE == 0
  #define INA239_CONFIG_SET_ADCRANGE  (0x00 << 4)  // Set bit 4 LOW in Register CONFIG
#elif INA239_ADC_RANGE == 1
  #define INA239_CONFIG_SET_ADCRANGE  (0x01 << 4)  // Set bit 4 HIGH in Register CONFIG
#endif
                                                  
// ADC_CONFIG sets the conversion times of the voltage, shunt voltage and
// temperature measurements as well as the ADC sample averaging count 
// see datasheet pages 21 and 29f.
// measurements to 50 µs; ADC sample averaging count to 1
#define INA239_ADC_CONFIG_FAST  (0x0F << 12 | 0x00 << 9 | 0x00 << 6 | 0x00 << 3 | 0x00)
// measurements to 4120 µs; ADC sample averaging count to 1024
#define INA239_ADC_CONFIG_SLOW  (0x0F << 12 | 0x07 << 9 | 0x07 << 6 | 0x07 << 3 | 0x07)
// measurements to 1052 µs; ADC sample averaging count to 4
#define INA239_ADC_CONFIG_SWEET (0x0F << 12 | 0x05 << 9 | 0x05 << 6 | 0x05 << 3 | 0x01)


class INA239 {
  public:
    
    // Resets all registers to default values. Setting the 'Reset Bit' to '1' 
    // generates a system reset that is the same as power-on reset.
    void reset() {
      RegisterWrite(INA239_WRITE_CONFIG, INA239_CONFIG_RESET);
    }

    void setup() {
      RegisterWrite(INA239_WRITE_CONFIG, INA239_CONFIG_SET_ADCRANGE);
      RegisterWrite(INA239_WRITE_ADC_CONFIG, INA239_ADC_CONFIG_SWEET);
      RegisterWrite(INA239_WRITE_SHUNT_CAL, INA239_SHUNT_CAL);
      delayMicroseconds(300);                                                   // TODO: is this too much? Look at datasheet how long

      #if DEBUG    
        int16_t config, config_adc, shunt;
        
        config = RegisterRead(INA239_READ_CONFIG);          // Relevant data in bits 15-0
        config_adc = RegisterRead(INA239_READ_ADC_CONFIG);  // Relevant data in bits 15-0
        shunt = RegisterRead(INA239_READ_SHUNT_CAL);        // Relevant data in bits 14-0
        
        Serial.print("CONFIG: ");
        Serial.print(config, BIN);
        Serial.print(", ");
        Serial.println(config, HEX);
        Serial.print("ADC_CONFIG: ");
        Serial.print(config_adc, BIN);
        Serial.print(", ");
        Serial.println(config_adc, HEX);
        Serial.print("SHUNT_CAL: ");
        Serial.print(shunt, BIN);
        Serial.print(", ");
        Serial.println(shunt, HEX);
      #endif
    }

    // Differential voltage in mV measured across the shunt output.
    float get_shunt_voltage() {
      int16_t data;
      float voltage;

      data = RegisterRead(INA239_READ_VSHUNT);  // Relevant data in bits 15-0
      voltage = data * INA239_VSHUNT_LSB;       // Conversion factor

      #if DEBUG    
        Serial.println("Try to get voltage");
        Serial.println(data);
        Serial.print(voltage);
        Serial.println(" mV");
      #endif

      return voltage;
    }

    // Voltage at PIN8, in V. Typically, the supply voltage of the MC.
    float get_pin8_voltage() {
      int16_t data;
      float voltage;

      data = RegisterRead(INA239_READ_VOLTAGE);  // Relevant data in bits 15-0
      voltage = data * INA239_Voltage_LSB;       // Conversion factor

      #if DEBUG    
        Serial.println("Try to get voltage");
        Serial.println(data);
        Serial.print(voltage);
        Serial.println(" V");
      #endif

      return voltage;
    }

    // Internal die temperature measurement. Two's complement value.
    // Temperature from -40 °C to +125 °C.
    float get_temperature() {
      int16_t data;
      float temperature;

      data = RegisterRead(INA239_READ_TEMPERATURE) >> 4;  // Relevant data in bits 15-4
      temperature = data * INA239_Temperature_LSB;        // Conversion factor

      #if DEBUG    
        Serial.println("Try to get temperature");
        Serial.println(data);
        Serial.print(temperature);
        Serial.println(" °C");
      #endif

      return temperature;
    }

    // Calculated current output in Amperes. Two's complement value.
    float get_current() {
      int16_t data;
      float current;

      data = RegisterRead(INA239_READ_CURRENT);  // Relevant data in bits 15-0
      current = data * INA239_Current_LSB;       // Conversion factor

      #if DEBUG    
        Serial.println("Try to get current");
        Serial.println(data);
        Serial.print(current);
        Serial.println(" A");
      #endif

      return current;
    }

  private:
    // Start an I/O operation for the INA239
    void startIO() {
      SPI.beginTransaction(INA239_SPI_SETTINGS);
      digitalWrite(INA_1_SPI_CS, LOW);
    }

    // Finish an I/O operation for the INA239
    void endIO() {
      digitalWrite(INA_1_SPI_CS, HIGH);
      SPI.endTransaction();
    }

    // Read specific INA239 Register
    int16_t RegisterRead(byte frame) {
      int16_t data;

      startIO();
      SPI.transfer(frame);
      data = SPI.transfer16(0x0);
      endIO();

      #if DEBUG    
        Serial.print("Try to read from Register ");
        Serial.print(frame >> 2, HEX);
        Serial.print("h => SPI Frame: ");
        Serial.println(frame, BIN);
        Serial.print("Data: ");
        Serial.println(data, BIN);
      #endif
      
      return data;
    }

    // Write 'command' to specific INA239 Register
    void RegisterWrite(byte frame, uint16_t command) {
      startIO();
      SPI.transfer(frame);
      SPI.transfer16(command);
      endIO();

      #if DEBUG    
        Serial.print("Try to write to Register ");
        Serial.print(frame >> 2, HEX);
        Serial.print("h => SPI Frame: ");
        Serial.println(frame, BIN);
        Serial.print("command: ");
        Serial.println(command, BIN);
      #endif
    }

};

// Create two global instances of class 'INA239'
INA239 INA_1;
INA239 INA_2;

#endif    // _VIPERLEED_INA239