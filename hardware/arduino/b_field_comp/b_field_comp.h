/*
ViPErLEED - Magnetic-field compensation
---------------------
Authors: Michele Riva, Christoph Pfungen, Stefan Mitterhöfer
16.05.2023, MR, CP: First version
                        Definition of Arduino Pins
                        Class MotorDriver
                        CLass Coil
27.11.2024, SM: Overhaul whole code
                        Move Class MotorDriver to separate file
                        Wrap INA239 and TimerCounter methods into Coil
---------------------
*/

#ifndef _VIPERLEED_B_FIELD_COMP
#define _VIPERLEED_B_FIELD_COMP


/** ------------------------- Communication with PC ------------------------ **/

// PC requested to set a specific duty cycle (ASCII 'D')
#define PC_SET_DUTY_CYCLE         68
// PC requested to set a specific current (ASCII 'I')
#define PC_SET_CURRENT               73
// PC requested measurement (ASCII 'M')
#define PC_MEASURE                  77
// PC requested time constant (ASCII 'T')
#define PC_SET_TIME_CONSTANT          84
// PC requested enable dynamic (ASCII 'E')
#define PC_SET_ENABLE_DYNAMIC          69
// PC requested to set transformation matrix (ASCII 't')
#define PC_SET_TRANSFORMATION_MATRIX          116
// PC requested to set calibration curve (ASCII 'c')
#define PC_SET_CALIBRATION_CURVE          99
// PC_SET_UP_ADCS: PC requested to prepare the
// ADCs for a measurement (ASCII 'S').
#define PC_SET_UP_ADCS       83

/** ------------------------- Finite state machine ------------------------- **/
// Set duty cycle for the coil using PWM.
#define STATE_SET_DUTY_CYCLE      50
// Set current for the coil.
#define STATE_SET_CURRENT          51
// Measure the current and voltage of the respective coil.
#define STATE_MEASURE              52
// Set the time constant.
#define STATE_SET_TIME_CONSTANT          53
// Enable/disable dynamic magnetic field compensation.
#define STATE_SET_ENABLE_DYNAMIC      54
// Set transformation matrix.
#define STATE_SET_TRANSFORMATION_MATRIX      55
// Set calibration curve.
#define STATE_SET_CALIBRATION_CURVE  56
// STATE_SET_UP_ADCS: Pick correct ADC channels and no. of measurement points.
#define STATE_SET_UP_ADCS          1


// The pins belonging to 'COIL_1_PWM' and 'COIL_2_PWM' should not be changed.
// Changing these may require choosing a different Timer/Counter module.
// Possible Pins for TimerCounter1: 9, 10, 11
// Possible Pins for TimerCounter4: 6, 10, 13
#define COIL_1_PWM              10   // PWM output 1, i.e., voltage value; Also A10; PB6 on Atmega32U4
#define COIL_1_ENABLE           2    // Used for motor driver; PD1 on ATmega32U4

#define COIL_2_PWM              6   // PWM output 2, i.e., voltage value; Also A7; PD7 on Atmega32U4
#define COIL_2_ENABLE           3    // Used for motor driver; PD0 on ATmega32U4

// Current direction: positive or negative?
#define COIL_1_SIGN             4    // Used for motor driver; PD4 on ATmega32U4
#define COIL_2_SIGN             5    // Used for motor driver; PC6 on ATmega32U4

#define INA_1_SPI_CS            7   // Pin for SPI communication with measurement device; PE6 on ATmega32U4
#define INA_2_SPI_CS            8   // Pin for SPI communication with measurement device; PB4 on ATmega32U4


// Generic error codes
enum error_t {
    NoError,
    OutOfRange,
    InvalidIOPin,
    InvalidChannel,
    InvalidPrescaler,
    NotImplemented,
};

#include "INA239.h"         // For configuration and functions of current measurement device INA239
#include "motordriver.h"    // for motor driver functions
#include "pwm.h"            // for set_pwm_frequency, set_coil_current

#if TLE7209_USE_SPI
  #define COIL_1_SPI_CS           11   // For SPI communication with the motordriver; PB7 on ATmega32U4
  #define COIL_2_SPI_CS           12   // For SPI communication with the motordriver; PD6 on ATmega32U4
#else
  #define COIL_Error              9   // MotorDriver Error-Flag (Error on LOW)
#endif


// Include some custom libraries. They are in the ../lib directory.
// This does NOT WORK CORRECTLY when using the Arduino IDE, as the
// latter copies all files into a temporary directory before compiling.
// If you want to use the Arduino IDE, you will have to manually
// Sketch>Add file... and include the following. Bear in mind, however,
// that this DUPLICATES the files in the b_field_comp folder. Thus any
// edit in the original library files WILL NOT BE REFLECTED in the files
// added manually. As an alternative, use upload_sketch.py for compiling,
// and/or uploading, as that one automatically takes care of linking
// library root directories before compilation.


class Coil {
    public:
        INA239 measurement;
        TimerCounter *pwm;
        MOTORDRIVER driver;

        // Class constructor, including member initializer list
        Coil(TimerCounter *pwm, byte enable, byte spi_cs_ina)
        : driver(MOTORDRIVER(enable)),
          measurement(INA239(spi_cs_ina)),
          pwm(pwm) {}

        // Class constructor, including member initializer list
        Coil(TimerCounter *pwm, byte spi_cs_motordriver, byte enable, byte spi_cs_ina)
        : driver(MOTORDRIVER(spi_cs_motordriver, enable)),
          measurement(INA239(spi_cs_ina)),
          pwm(pwm) {}

        void setup() {
          #if DEBUG
            Serial.println("Coil setup started ...");
          #endif

          measurement.reset();
          measurement.setup();
          pwm->setup();
          driver.reset();
          driver.setup();

          #if DEBUG
            Serial.println("Coil setup finished.");
          #endif                                                                    // TODO: Implement Self-Calibration and Degausing Routine in setup()
        }


//        byte set_current(float coil_current) {                                 // TODO: Implement set_current()
//            byte err = set_signed_pwm_value(coil_current, pwm_sign_pin,
//                                            tc4_reg_addr);
//            if(err)
//                return err;
//            _current_setpoint = coil_current;
//        }

        // just here for tests, later in sec. private
        error_t set_duty_cycle(float value) {
            /**Set PWM duty cycle (average voltage), including sign output.

            Parameters
            ----------
            value : float
                PWM duty cycle to be set. Should be between -1.0 and +1.0.

            Returns
            -------
            error_code : byte
                0 for no error
                1 for 'value' (aka coil current) out-of-range
            **/
            error_t err = pwm->set_duty_cycle(value);
            if(err)
                return err;

            return NoError;
        }
        // -------------------

        float get_current() {
        // Returns float: current in Amperes
            return measurement.get_current();
        }

        float get_voltage() {
        // Returns float: voltage in Volt
            return measurement.get_pin8_voltage();
        }

        float get_ambient_temperature() {
        // Returns float: temperature in °C
            return measurement.get_temperature();
        }
                                                                                // TODO: Methode for coil: adjust_current_to_setpoint(), which is called in main-loop

        error_t get_characteristic_curve(float delta) {
            // Get the characteristic curves (current/voltage over duty cycle)
            // of the system (Arduino + MD + Coil) and send it via serial
            // interface to the host computer.
            //
            // Parameters
            // ----------
            // delta : float
            //   The difference between two setpoints of the duty cycle
            //
            // Returns
            // -------
            // error_code : byte
            //    0 for no error
            //    1 for 'value' (aka coil current) out-of-range
            // -------

            error_t err = NoError;
            float voltage, current;
            float setpoint = -1.000;

            Serial.setTimeout(100);
            Serial.begin(9600);  // Open serial port, set data rate to 9600 bps
            delay(2000);         // The IDE needs some time to connect to the Serial Monitor

            Serial.println("duty cycle;voltage;current");
            do {
                err = set_duty_cycle(setpoint);
                if(err)
                    return err;
                delay(100);
                voltage = get_voltage();
                current = get_current();

                Serial.print(setpoint,3);
                Serial.print(";");
                Serial.print(voltage,3);
                Serial.print(";");
                Serial.println(current,6);

                setpoint += delta;
            } while(abs(1 - setpoint) > 1E-10);

            return NoError;
        }

        void set_calibration_values() {


        }

    private:
        float _current_setpoint;
};

// Create two global instances of class 'TimerCounter' with respective initializers
TimerCounter1 pwm1(COIL_1_PWM, COIL_1_SIGN);
TimerCounter4 pwm2(COIL_2_PWM, COIL_2_SIGN);
// Create two global instances of class 'Coil' with respective initializers
#if TLE7209_USE_SPI
  Coil coil_1(&pwm1, COIL_1_SPI_CS, COIL_1_ENABLE, INA_1_SPI_CS);
  Coil coil_2(&pwm2, COIL_2_SPI_CS, COIL_2_ENABLE, INA_2_SPI_CS);
#else
  Coil coil_1(&pwm1, COIL_1_ENABLE, INA_1_SPI_CS);
  Coil coil_2(&pwm2, COIL_2_ENABLE, INA_2_SPI_CS);
#endif

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

// Two floats that represent the duty cycles of each respective coil.
floatOrBytes dutyCycle[2];
// Two current values that have to be converted to duty cycles.
floatOrBytes targetCurrent[2];
// Two floats which represent the time constants of the coils.
floatOrBytes timeConstant[2];
// Enable/disable dynamic magnetic field compensation.
bool enableDynamic;
// Transformation matrix which converts external field to field on sample.
floatOrBytes transformationMatrix[3][2];
// Calibration curve for converting PWM duty cycles to currents.
floatOrBytes calibrationCurve[2][12];
// Resistance of the coils in ohms.
float coilResistance[2];

/**
bool measureResistance1 = false;
bool measureResistance2 = false;
**/
float current_max;
unsigned long current_time = 0;
unsigned long last_current_time = 0;

float new_current[2];
float last_current[2];
float new_voltage[2];

bool measured = false;

#endif  // _VIPERLEED_B_FIELD_COMP
