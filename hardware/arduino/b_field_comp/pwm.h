/*
ViPErLEED - Driver for using fast-PWM mode via Timer/Counter.
---------------------
Authors: Michele Riva, Christoph Pfungen, Stefan Mitterhöfer
15.05.2023, MR, CP: First version with TC4
17.07.2024, SM: Overhaul whole code as class, include TC1, change introduction
---------------------
*/

/**
Quick Start Guide to Timer/Counter
==================================

The ATmega32U4 includes the following Timer/Counter modules:
- 8-bit Timer/Counter0
- 16-bit Timers/Counters (TC1 and TC3)
- 10-bit High Speed Timer/Counter4

This section is intended as a concise introduction to the PWM waveform
generation using the Timer/Counter1 or Timer/Counter4, subsequently abbreviated
as TCn, where n = 1 or 4.

The TCn counts up from BOTTOM to TOP. When reaching TOP TCn either restarts
at BOTTOM or counts down to BOTTOM again (depending on the configuration).
The difference of BOTTOM and TOP determines the frequency of the PWM signal.
TCn contains independant Output Compare Units OCnx, 
where x = A, B or C for TC1 and x = A, B, C or D for TC4 respectively.
While counting up (or down) these Compare Units OCnx are constantly comparing
their register values OCRnx to the TCn register values TCNTn. If a match occurs
defined Output Pins are either set or cleared (depending on the configuration)
to generate the PWM signal. Only a few Pins on the Arduino Board can be 
connected to the OCnx (see datasheet Figure 1-1 and section 10.3).
The value of the Compare Unit OCnx determines the Duty cylce of the PWM signal.

The relevant formulae are also given in the datasheet sec. 14.8.3 and 15.8.
For TCn operated in Fast PWM Mode the PWM frequency, TCn clock frequency and 
PWM resolution are determined as follows:

(1)   f_PWM = f_clk_TCn / (TOP + 1)

(2)   f_clk_TCn = F_CPU_CLK / TCn_prescaler

(3)   res_PWM = log2(TOP + 1) = log2(f_clk_TCn / f_PWM)

Every time TCNTn has a match with TOP, TCNTn is reset to zero (or set to 1
respectively) and begins to count up or down again until another match occurs.
The TCn frequency 'f_clk_TCn' is derived from the CPU clock 'F_CPU_CLK' 
divided by the TCn prescaler settings. If a value smaller than three is written
to TOP, the hardware will replace that value with three in order to prevent
too high PWM frequency.

Equation (3) shows that the user may want to stay as close to the maximum value
of TOP as possible, in order to achieve the highest possible PWM resolution.
Therefore, 'set_pwm_frequency' will first determine the TCn prescaler value
which will allow for the highest possible 10/16-bit value inside TOP which
in turn will provide the best possible PWM resolution. 


---------- Timer/Counter1 ----------
To use Fast PWM and the possibility to change the value of TOP while operation
'Mode' is set to 15. Therefore the bits WGM13:WGM10 are set to 1 (via Registers
TCCR1A and TCCR1B). Detailed description see datasheet section 14.10 and 
table 14-4. As the TOP value is stored in Register OCR1A only OC1B and OC1C are
available to generate the PWM signal. Thus either Pin 10/PB6 (OC1B) or 
Pin 11/PB7 (OC1C) can be used on the Arduino Board to generate a PWM signal
using TC1 (see datasheet Figure 1-1 and section 10.3.1). To activate Pin 10/PB6
the corresponding Bit 'DDB6' in Register DDRB has to be set to 1 (see datasheet
sections 10.3.1 and 10.4.3). Analogous for Pin 11/PB7 (Bit 'DDB7' in Register
DDRB to 1).

An inverted or non-inverted PWM signal can be archieved by setting the 
'Compare Output Mode' bits COM1x1:COM1x0 (via Register TCCR1A) to 0b10 
(clear at match, set at TOP) or 0b11 (set at match, clear at TOP), respectively.

Prescaling is adjusted by the Clock Select bits CS12:CS10 (via Register TCCR1B).
Details see datasheet table 14-5. This feature is not used for TC1, because a
frequency less than 244 Hz (16 MHz/2^16) is usually not needed. So the
prescaler '1' is initialized during setup.


---------- Timer/Counter4 ----------
In contrast to the other TC modules, TC4 can be operated with up to 64 MHz if
the supply voltage is higher than 4V. The firmware does not make use of this 
feature.

Note: The 11-bit Enhanced PWM Mode is not working up to some chip revision.

To use Fast PWM and the possibility to change the value of TOP while operation
'Waveform Generation Mode' is set to 0. Therefore the bits WGM41:WGM40 are set
to 0b00 (via Register TCCR4D) and bit PWM4X is set to 1 (via Register TCCR4B).
Detailed description see datasheet section 15.12 and table 15-18.
As the TOP value is stored in Register OCR4C only OC4A, OC4B and OC4D are
available to generate a PWM signal. Thus either Pin 13/PC7 (OC4A),
Pin 10/PB6 (OC4B) or Pin 6/PD7 (OC4D) can be used on the Arduino Board to
generate a PWM signal using TC4 (see datasheet Figure 1-1 and sections 10.3.1,
10.3.2 and 10.3.3). To activate Pin 13/PC7 the corresponding Bit 'DDC7' in
Register DDRC has to be set to 1 (see datasheet sections 10.3.2 and 10.4.6).
Analogous for Pin 10/PB6 (Bit 'DDB6' in Register DDRB to 1, see datasheet
sections 10.3.1 and 10.4.3) and Pin 6/PD7 (Bit 'DDD7 in Register DDRD to 1,
see datasheet sections 10.3.3 and 10.4.9).

Depending on the 'Compare Output Mode' settings COM4x1:0 (via Registers TCCR4A
and TCCR4C), the Output Compare pins OC4x will be set or cleared differently.
Set COM4x1:0 to 0b10 to clear at match and set when TCNT4 = 0x000
(0b11: set at match, clear when 0x000).
                                                   
Prescaling is adjusted by the Clock Select bits CS43:CS40 (via Register TCCR4B).
Details see datasheet table 15-14.

TC4 can be used in 8-bit or 10-bit operation, depending on whether the user 
also programs the two high bits inside shared register TC4H (cf. datasheet 
sec. 15.2.2). For example, in order to perform the 10-bit write 0x31F to 
register TC4H:OCR4C, first program TC4H with 0x3, then write 0x1F to OCR4C.


**/

#ifndef _VIPERLEED_B_FIELD_PWM
#define _VIPERLEED_B_FIELD_PWM

#define F_CPU_CLK 16e6                  // Arduino Micro CPU clock = 16 MHz

// Min. and max. PWM frequency see equation (1) and (2)
#define TC1_F_PWM_MIN (F_CPU_CLK / (1024 * (pow(2, 16) + 1)))  // F_PWM_MIN = 0,238 Hz
#define TC1_F_PWM_MAX (F_CPU_CLK / (1 * (3 + 1)))              // F_PWM_MAX = 4 MHz

#define TC4_F_PWM_MIN (F_CPU_CLK / (16384 * (pow(2, 10) + 1))) // F_PWM_MIN = 0,953 Hz
#define TC4_F_PWM_MAX (F_CPU_CLK / (1 * (3 + 1)))              // F_PWM_MAX = 4 MHz

// Addresses of TCn registers
#define TC1_ADDR        &TCNT1
#define TC1_ADDR_TOP    &OCR1A
#define TC1_MAX         0xFFFF          // see datasheet section 14.1.2

#define TC4_ADDR        &TCNT4
#define TC4_ADDR_TOP    &OCR4C
#define TC4_MAX         0x3FF           // see datasheet section 15.2.5


class TimerCounter {        
    public:
        TimerCounter(byte pwm, byte sign) 
        : _PWM_PIN(pwm), _PWM_SIGN_PIN(sign) {}

        virtual void setup() = 0;
        virtual void tc_sfr_reset() = 0;
//        Put here all methods, which should be available for class 'Coil' 
        
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
            
            Notes
            -----
            First change the PWM frequency, if necessary.
                For duty cycle values, which will generate shorter PWM pulses 
                than 'MD_PWM_PULSE_MIN', the PWM frequency will be decreased 
                automatically to avoid to high discrepancy in the set and
                actual value of the current. See 'set_pwm_frequency'
            Then set the duty cycle.
                Notice that the 'coil_current', i.e. the time-averaged value 
                of the signal from the PWM, is exactly the same as the duty 
                cycle of the PWM itself. Set the PWM duty cycle at 
                '_tc_reg_addr': The duty cycle is the desired 'on' time in 
                percent, i.e. 'coil_current' times the PWM period measured in 
                TCn clock ticks.
            Finally set the sign depending on the current direction.
            **/
            
            float abs_value;
            float pwm_frequency;
            bool change_sign = false;
            
            // If duty cycle hasn't changed, don't do anything.
            if (abs(value - _duty_cylce) < 1E-10) {
                return NoError;
            }
            
            if (value < -1.0 || value > 1.0) 
                return OutOfRange;
            
            // Prepare change of direction, 
            // if it has changed compared to last value
            if (sgn(_duty_cylce) != sgn(value)) {                           
                change_sign = true;                                                                             
            }
            // As sign is already prepared, just use absolute value from here.
            abs_value = abs(value);
            
            #if DEBUG
              Serial.print("Last Duty Cylce _duty_cylce = ");
              Serial.println(_duty_cylce);
              Serial.print("New Duty Cycle = ");
              Serial.println(value);
            #endif
            _duty_cylce = value;
            
            if (abs_value > 1E-10) {
                // Change PWM frequency, if necessary   
                pwm_frequency = find_optimum_frequency(abs_value);
                set_pwm_frequency(pwm_frequency);                                       
            }

            set_duty_cycle_to_register(round(abs_value * (_pwm_period)));

            // change of direction, if it has changed compared to last value
            if (change_sign) {                                                  // TODO: Disable MotorDriver while switching from neg. to pos. or pos. to neg. ?
                set_sign(value);                                                //       (just if necessary, probably not)
            }
            
            return NoError;
        }

    protected:
        enum TC_PWM_CHANNEL { TC_PWM_CH_A, TC_PWM_CH_B,
                              TC_PWM_CH_C, TC_PWM_CH_D };
        const byte _PWM_PIN;
        const byte _PWM_SIGN_PIN;
        float _TC_F_PWM_MIN;
        float _TC_F_PWM_MAX;
        TC_PWM_CHANNEL _tc_pwm_channel; // Channel for PWM duty cycle
        // Period and frequency are stored to avoid divisions during runtime.
        uint16_t _pwm_period;           // currently PWM periode
        float _f_pwm = -1;              // currently PWM frequency
        float _duty_cylce = -2;         // currently PWM duty cycle
        
        virtual TC_PWM_CHANNEL pin_to_tc_channel() = 0;
        virtual error_t enable_pwm_channel(bool) = 0;
        virtual error_t connect_pwm_output(byte) = 0;
        virtual void enable_fast_pwm_mode() = 0;
        virtual error_t set_polarity_inverted(bool) = 0;
        virtual void set_pwm_period(float) = 0;
        virtual void set_duty_cycle_to_register(uint16_t) = 0;
//         virtual uint16_t find_optimum_prescaler(float) {}    // NOT USED
//         virtual error_t set_pwm_clock_prescaler(uint16_t) {} // NOT USED

        void set_value(uint16_t value, uint16_t *reg) {
            /**Write 16-bit value to Timer/Counter1 register.
        
            Parameters
            ----------
            value : uint16_t
                value to be written
            reg : uint16_t *
                Address of register to be written to
        
            Notes
            -----
            Registers are 8 bits. To write 16-bit values send the HIGH Byte
            first. Note that when using “C”, the compiler handles the 16-bit
            access. More info: Atmega32U4 datasheet, section 14.2.
            **/
            
            #if DEBUG
              Serial.print("Write value ");
              Serial.print(value);
              Serial.print(" (");
              Serial.print(value,BIN);
              Serial.print(") to 16-bit register ");
              Serial.println((long) reg);
            #endif
        
            *reg = value;
        }
        
        void set_value(uint16_t value, uint8_t *reg) {
            /**Write 10-bit value to Timer/Counter4 register.
        
            Parameters
            ----------
            value : uint16_t
                value to be written
            reg : uint8_t *
                Address of register to be written to
        
            Notes
            -----
            Registers are 8 bits. Some can also accept 10-bit values.
            To do this, the two highest bits are to be written in the
            *shared* register TC4H right before the remaining 8 bits
            are written in the desired register. More info: Atmega32U4
            datasheet, section 15.11.
            **/
            
            #if DEBUG
              Serial.print("Write value ");
              Serial.print(value);
              Serial.print(" (");
              Serial.print(value,BIN);
              Serial.print(") to 10-bit register ");
              Serial.println((long) reg);
            #endif
        
            noInterrupts();  // Ensure nothing bothers setting two registers
            TC4H = value >> 8;
            *reg = value & 255;
            interrupts();
        }
        
        void set_sign(float value) {
            /**Set digital I/O pin according to current direction in a coil.
            Also set the PMW polarity.
        
            Parameters
            ----------
            value : float
                Coil current direction: Either positive or negative
            **/

            if (value < 0) {                
                set_polarity_inverted(true);
                digitalWrite(_PWM_SIGN_PIN, HIGH);
            }
            else {
                set_polarity_inverted(false);
                digitalWrite(_PWM_SIGN_PIN, LOW);
            }
        }
        
        float find_optimum_frequency(float value) {
            /**Find the optimum PWM frequency.
            For duty cycle values, which will generate shorter PWM pulses 
            than 'MD_PWM_PULSE_MIN', the PWM frequency will be decreased 
            automatically to avoid to high discrepancy in the set and
            actual value of the current. See 'set_pwm_frequency'
            
            Parameters
            ----------
            value : float
                Duty cycle of the PWM signal.
        
            Returns
            -------
            optimum PWM frequency: float
            **/
            
            if (value < DUTY_CYCLE_MIN) {
                //pwm_pulse = value / F_PWM; // needed PWM pulse duration [s]
                //pwm_frequency = (F_PWM / MD_PWM_PULSE_MIN) * pwm_pulse;    
                //pwm_frequency = value / MD_PWM_PULSE_MIN;
                return value * MD_PWM_FREQUENCY_MAX;
            }    
            else if (value > DUTY_CYCLE_MAX) {
                //pwm_frequency = (1 - value) / MD_PWM_PULSE_MIN;
                return (1 - value) * MD_PWM_FREQUENCY_MAX;
            }                   
            else {
                return F_PWM;
            }
        }
        
        error_t set_pwm_frequency(float value) {
            /**Set PWM frequency as period in register TOP.
        
            Parameters
            ----------
            value : float
                PWM output frequency in Hz.
                Notice that the resolution of the PWM scales with PWM frequency.
                The PWM frequency is shared across all available 
                Timer/Counter channels A, B, C and D.
        
            Returns
            -------
            error_code : byte
                0 for no error
                1 for requested pwm frequency out-of-range
            **/
            
            // If PWM frequency hasn't changed, don't do anything.
            if (abs(value - _f_pwm) < 1E-10) {
                return NoError;
            }
            
            // If PWM frequency is out of range,
            // correct it to MIN or MAX, respectively.
            if (value < _TC_F_PWM_MIN) {
                value = _TC_F_PWM_MIN;
            }
            else if (value > _TC_F_PWM_MAX) {
                value = _TC_F_PWM_MAX;
            }
            
            #if DEBUG
              Serial.print("Last PWM frequency _f_pwm = ");
              Serial.println(_f_pwm);
              Serial.print("New PWM frequency value = ");
              Serial.println(value);
            #endif
            _f_pwm = value;
            
            set_pwm_period(value);
            
            return NoError;
        }
        
        error_t set_config() {
            /**Set Arduino configuration
                Enable selected pins and consequential channels
                Enable fast PWM mode
        
            Returns
            -------
            error_code : byte
                0 for no error
                2 for invalid I/O pin
                3 for invalid TC channel
            **/
            error_t err = NoError;
            
            _tc_pwm_channel = pin_to_tc_channel();
            
            #if DEBUG
              Serial.print("_PWM_PIN = ");
              Serial.println(_PWM_PIN);
              Serial.print("_PWM_SIGN_PIN = ");
              Serial.println(_PWM_SIGN_PIN);
              Serial.print("_tc_pwm_channel = ");
              Serial.println(_tc_pwm_channel);
            #endif        
    
            pinMode(_PWM_PIN, OUTPUT);
            pinMode(_PWM_SIGN_PIN, OUTPUT);
            
            enable_fast_pwm_mode();
            err = enable_pwm_channel(true);
            if(err)
                return err;
            err = connect_pwm_output(_PWM_PIN);
            
            return err;                                     
        }
        
        error_t PWM_init() {  
             /**Set PWM frequency and duty cycle to initial values.
             
            Returns
            -------
            error_code : byte
                0 for no error
                1 for out-of-range
            **/
            error_t err = NoError;
            
            err = set_pwm_frequency(F_PWM);
            if(err)
                return err;
            err = set_duty_cycle(0.0);

            return err;
        }
};

class TimerCounter1 : public TimerCounter {
    public:   
        TimerCounter1(byte pwm, byte sign) 
        : TimerCounter(pwm, sign) {}
        
        void setup() override {
            #if DEBUG
              Serial.println("TC1 setup started ...");
            #endif
            
            tc_sfr_reset();
            
            _TC_F_PWM_MIN = TC1_F_PWM_MIN;
            _TC_F_PWM_MAX = TC1_F_PWM_MAX;
            _TOP = TC1_ADDR_TOP;
            _TCNTn = TC1_ADDR;
            _MAX = TC1_MAX;
            
            _tc_reg_addr = pin_to_tc_reg_addr();
            TCCR1B |= 0b001<<CS10;       // set prescaler to 1          
            
            set_config();
            PWM_init();
            
            #if DEBUG
              Serial.print("_TC_F_PWM_MIN = ");
              Serial.println(_TC_F_PWM_MIN);
              Serial.print("_TC_F_PWM_MAX = ");
              Serial.println(_TC_F_PWM_MAX);
              Serial.print("_TOP = ");
              Serial.println(*_TOP);
              Serial.print("_TCNTn = ");
              Serial.println(*_TCNTn);
              Serial.print("_MAX = ");
              Serial.println(_MAX,HEX);
              Serial.print("_tc_reg_addr = ");
              Serial.println(*_tc_reg_addr);
              Serial.println("TC1 setup finished.");
            #endif
        }

        void tc_sfr_reset() override {
            // Some special function registers belonging to
            // Timer/Counter1 are not properly zero'ed on POR
            TCCR1A = 0;
            TCCR1B = 0;
            TCCR1C = 0;
            OCR1AH = 0;
            OCR1AL = 0;
            OCR1BH = 0;
            OCR1BL = 0;
            OCR1CH = 0;
            OCR1CL = 0;
            TIMSK1 = 0;
            TIFR1 = 0;
        }
        
    private:
        uint16_t *_tc_reg_addr;         // Register address for PWM duty cycle
        uint16_t *_TOP;                 // Register address for PWM periode
        uint16_t *_TCNTn;               // Register address of Counter n
        uint16_t _MAX;                  // Maximum value of Counter n
        
        unsigned int read_register_16bit(uint16_t *reg ) {
            unsigned int i;
            /* Read *reg into i */
            i = *reg;
            
            return i;
        }
        
        TC_PWM_CHANNEL pin_to_tc_channel() override {
            switch(_PWM_PIN) {
              case  9: return TC_PWM_CH_A;
              case 10: return TC_PWM_CH_B;
              case 11: return TC_PWM_CH_C;
              default: return 0;                                                       // Change to more descriptive error
            }
        }
        
        uint16_t* pin_to_tc_reg_addr() {
            switch(_PWM_PIN) {
              case  9: return _SFR_ADDR(OCR1A);
              case 10: return _SFR_ADDR(OCR1B);
              case 11: return _SFR_ADDR(OCR1C);
              default: return 0;                                                       // Change to more descriptive error
            }
        }
        
        error_t enable_pwm_channel(bool enable) override {
            /**Enable/disable PWM output on TC channels A, B or C.
        
            Parameters
            ----------
            enable : bool
                Enable or disable the channel
        
            Returns
            -------
            error_code : error_t
                0 for no error
                3 for invalid TC channel
        
            **/    
            switch (_tc_pwm_channel) {
              case TC_PWM_CH_A: 
                if(enable) {
                  TCCR1A |= (1<<COM1A1);
                  TCCR1A &= ~(1<<COM1A0);     // set bit 'COM1A0' to 0
                }
                else
                  // set bits 'COM1A1:0' to 0
                  TCCR1A &= ~((1<<COM1A1)|(1<<COM1A0));
                break;
                
              case TC_PWM_CH_B: 
                if(enable) {
                  TCCR1A |= (1<<COM1B1);
                  TCCR1A &= ~(1<<COM1B0);     // set bit 'COM1B0' to 0
                }
                else
                  // set bits 'COM1B1:0' to 0
                  TCCR1A &= ~((1<<COM1B1)|(1<<COM1B0));
                break;
        
              case TC_PWM_CH_C: 
                if(enable) {
                  TCCR1A |= (1<<COM1C1);
                  TCCR1A &= ~(1<<COM1C0);     // set bit 'COM1C0' to 0
                }
                else
                  // set bits 'COM1C1:0' to 0
                  TCCR1A &= ~((1<<COM1C1)|(1<<COM1C0)); 
                break;
              default: return InvalidChannel;
            }
            
            #if DEBUG
                if(enable)
                  Serial.print("Enable");
                else
                  Serial.print("Disable");
                Serial.print(" PWM output on TC channel ");
                Serial.println(_tc_pwm_channel);
            #endif
            
            return NoError;
        }
        
        error_t connect_pwm_output(byte io_pin) override {
            /**Connect internal Waveform Output 'OCnx' to physical 'io_pin'.
        
            Parameters
            ----------
            io_pin : byte
                Physical pin which will output PWM-generated waveform 'OCnx'
        
            Returns
            -------
            error_code : error_t
                0 for no error
                2 for invalid I/O pin
        
            **/
            
            #if DEBUG
              Serial.print("Connect PWM output pin ");
              Serial.println(io_pin);
            #endif
            
            switch(io_pin) {
              case  9: DDRB |= (1<<DDB5); break;
              case 10: DDRB |= (1<<DDB6); break;
              case 11: DDRB |= (1<<DDB7); break;
              default: return InvalidIOPin;
            }    
            return NoError;
        }
        
        void enable_fast_pwm_mode() override {
            /**Enable Fast PWM Mode on TC1.
            TOP Counter stored in OCR1A, see Table 14-4.
            **/
        
            // Enable Fast PWM Mode (WGM13:10 = 0b1111)
            TCCR1A |= (1<<WGM11)|(1<<WGM10);
            TCCR1B |= (1<<WGM13)|(1<<WGM12);
        }         
        
        error_t set_polarity_inverted(bool invert) override {
            /**Set the polarity of the PWM signal.
            Parameters
            ----------
            invert : bool
                Enable or disable inverted polarity
                
            Returns
            -------
            error_code : error_t
                0 for no error
                3 for invalid TC channel
            
            Notes
            -----        
            Depending on the 'Compare Output Mode' settings COM1x1:0 (via
            Register TCCR1A) the Output Compare pins OC1x will be set or
            cleared differently. Set COM1x1:0 to 0b10 to clear at match and
            set at TOP (0b11: set at match, clear at TOP).
            To avoid a voltage peak while changing the direction, due to
            change of _PWM_SIGN_PIN and PWM polarity are not in phase,
            set the counter TCNT1 to MAX value to reset it with the next tik.
            **/
            
            switch (_tc_pwm_channel) {
                case TC_PWM_CH_A:
                    if(invert)
                        TCCR1A |= (1<<COM1A1)|(1<<COM1A0);
                    else {
                        TCCR1A |= (1<<COM1A1);
                        TCCR1A &= ~(1<<COM1A0);     // set bit 'COM1A0' to 0
                    } 
                    break;
                
                case TC_PWM_CH_B:
                    if(invert) {
                        TCCR1A |= (1<<COM1B1)|(1<<COM1B0);
                    }
                    else {
                        TCCR1A |= (1<<COM1B1);
                        TCCR1A &= ~(1<<COM1B0);     // set bit 'COM1B0' to 0
                    }
                    break;
        
                case TC_PWM_CH_C:
                    if(invert)
                        TCCR1A |= (1<<COM1C1)|(1<<COM1C0);
                    else {
                        TCCR1A |= (1<<COM1C1);
                        TCCR1A &= ~(1<<COM1C0);     // set bit 'COM1C0' to 0
                    }    
                    break;
                default: return InvalidChannel;
            }
            
            set_value(_MAX, _TCNTn);
            
            #if DEBUG
                if(invert)
                    Serial.print("Inverted ");
                else
                    Serial.print("Non-inverted ");
              Serial.print("polarity on TC channel ");
              Serial.println(_tc_pwm_channel);
            #endif
            
            return NoError;                                                               
        }
        
        void set_pwm_period(float frequency) {
            /**Set PWM period and PWM clock prescaler.
            
            Parameters
            ----------
            frequency : float
                PWM output frequency in Hz.
                
            Notes
            -----
            Notice that the resolution of the PWM scales with PWM frequency.
            Using prescalers for lower frequencies is recommended.
            For 16-bit timer no prescaler is used, because only at frequencies 
            below 244 Hz (see eq. (3): 16 MHz / 2^16) a prescaler would be 
            necessary. As this frequencies are hardly ever used it is
            neglectable.
            **/
            
            // see eq. (1) and (2) with no prescaler (clk_ps = 1)
            _pwm_period = (F_CPU_CLK / frequency) - 1;    // TOP = _pwm_period
            
            set_value(_pwm_period, _TOP);     // TOP = _pwm_period
        }
        
        void set_duty_cycle_to_register(uint16_t value) override {
            set_value(value, _tc_reg_addr);
            
            #if DEBUG
              Serial.print("Set value ");
              Serial.print(value);
              Serial.print(" to register ");
              Serial.println(*_tc_reg_addr);
            #endif
        }
        

// -------------------------------- NOT USED --------------------------------

//         const uint16_t _TC_CLK_PRESCALER[5] = {1, 8, 64, 256, 1024};

//         uint16_t find_optimum_prescaler(float value) {
//             /**Find optimum clock prescaler.
//         
//             Parameters
//             ----------
//             value : float
//                 For calculating optimum clock prescaler
//         
//             Returns
//             -------
//             uint16_t : clk_ps
//                 Returns best prescaler
//         
//             Notes
//             -----   
//             Calculate the TC1 clock prescaler: Use the smallest possible
//             prescaler which will result in the highest achievable PWM
//             resolution, i.e. make OCR1A as large as possible.
//             For a given CPU and PWM frequency, choose a clock prescaler that
//             keeps the PWM resolution as close as possible to the maximum
//             of 16 bits.
//             For example, a 16 MHz CPU clock and 20 kHz PWM yields 800 steps
//             in the interval [0,799]. This corresponds to a resolution of
///            9.64 bits.
//             **/
//             
//             float pwm_resolution;
//             float f_clk_tc;
//             uint16_t clk_ps;
// 
// 
//             // Iterate through the available TC1 prescaler values:
//             for(uint8_t i = 0; i < 4; i++) {
//                      
//                 clk_ps = _TC_CLK_PRESCALER[i];                
//                 f_clk_tc = F_CPU_CLK / clk_ps;              // eq. (2)
//                 pwm_resolution = log2(f_clk_tc / value);    // eq. (3)
//                 
//                 // Jump out of the loop, if best prescaler value is found.
//                 if(pwm_resolution <= 16.0) break;
//             }
//             
//             #if DEBUG
//               Serial.print("Found ");
//               Serial.print(clk_ps);
//               Serial.print(" as optimum prescaler for pwm frequency ");
//               Serial.println(value);
//             #endif
//             
//             return clk_ps;
//         }
// 
//         error_t set_pwm_clock_prescaler(uint16_t tc_clock_prescaler) {
//             /**Set PWM clock prescaler.
//         
//             Parameters
//             ----------
//             tc_clock_prescaler : uint16_t
//                 Sets the requested TC clock prescaler
//         
//             Returns
//             -------
//             error_code : error_t
//                 0 for no error
//                 4 for invalid clock prescaler
//             
//             **/
//             
//             byte ps_select;
//         
//             // Convert the TC1 clock prescaler value to the respective 
//             // prescaler select entry, see ATmega32U4 datasheet table 14-5
//             switch(tc_clock_prescaler) {
//               case 1: ps_select = 1; break;
//               case 8: ps_select = 2; break;
//               case 64: ps_select = 3; break;
//               case 256: ps_select = 4; break;
//               case 1024: ps_select = 5; break; 
//               default: return InvalidPrescaler;
//             }
//             
//             // To reset the prescaler (or counter TCNTn) set pins
//             // CS12:CS10 to zero (No clock source, Timer/Counter stopped)
//             TCCR1B &= ~(0b111<<CS10);        // set bits 'CS12:CS10' to 0           
//             TCCR1B |= ps_select<<CS10;       // set prescaler
//             
//             #if DEBUG
//               Serial.print("Write ");
//               Serial.print(ps_select,BIN);
//               Serial.print(" to register TCCR1B for prescaler ");
//               Serial.println(tc_clock_prescaler);
//             #endif
//             
//             return NoError;          
//         }

// --------------------------------------------------------------------------               
};


class TimerCounter4 : public TimerCounter {
    public:        
        TimerCounter4(byte pwm, byte sign) 
        : TimerCounter(pwm, sign) {}
        
        void setup() override {
            #if DEBUG
              Serial.println("TC4 setup started ...");
            #endif

            tc_sfr_reset();
            
            _TC_F_PWM_MIN = TC4_F_PWM_MIN;
            _TC_F_PWM_MAX = TC4_F_PWM_MAX;
            _TOP = (uint8_t*)TC4_ADDR_TOP;
            _TCNTn = (uint8_t*)TC4_ADDR;
            _MAX = TC4_MAX;
            
            _tc_reg_addr = (uint8_t*)pin_to_tc_reg_addr();
            
            set_config();
            PWM_init();
            
            #if DEBUG
              Serial.print("_TC_F_PWM_MIN = ");
              Serial.println(_TC_F_PWM_MIN);
              Serial.print("_TC_F_PWM_MAX = ");
              Serial.println(_TC_F_PWM_MAX);
              Serial.print("_TOP = ");
              Serial.println(*_TOP);
              Serial.print("_TCNTn = ");
              Serial.println(*_TCNTn);
              Serial.print("_MAX = ");
              Serial.println(_MAX,HEX);
              Serial.print("_tc_reg_addr = ");
              Serial.println(*_tc_reg_addr);
              Serial.println("TC4 setup finished.");
            #endif
        }
        
        void tc_sfr_reset() override {
            // Some special function registers belonging to
            // Timer/Counter4 are not properly zero'ed on POR
            TCCR4A = 0;
            TCCR4B = 0;
            TCCR4C = 0;
            TCCR4D = 0;
            TCCR4E = 0;
            OCR4A = 0;
            OCR4B = 0;
            OCR4C = 0;
            OCR4D = 0;
            TIMSK4 = 0;
            DT4 = 0;
        }
        
    private:
        uint8_t *_tc_reg_addr;          // Register address for PWM duty cycle
        uint8_t *_TOP;                  // Register address for PWM periode
        uint8_t *_TCNTn;                // Register address of Counter n
        uint8_t _MAX;                   // Maximum value of Counter n
        
        unsigned int read_register_10bit(uint8_t *reg ) {
            unsigned int i;
            /* Read *reg into i */
            i = *reg;
            i |= ((unsigned int)TC4H << 8);
            
            return i;
        }
        
        TC_PWM_CHANNEL pin_to_tc_channel() override {
            switch(_PWM_PIN) {
              case 13: return TC_PWM_CH_A;
              case 10: return TC_PWM_CH_B;
              case  6: return TC_PWM_CH_D;
              default: return 0;                                                       // Change to more descriptive error
            }
        }
        
        uint8_t* pin_to_tc_reg_addr() {
            switch(_PWM_PIN) {
              case 13: return _SFR_ADDR(OCR4A);
              case 10: return _SFR_ADDR(OCR4B);
              case  6: return _SFR_ADDR(OCR4D);
              default: return 0;                                                       // Change to more descriptive error
            }
        }
        
        error_t enable_pwm_channel(bool enable) override {
            /**Enable/disable PWM output on TC channels A, B or D.

            Parameters
            ----------
            enable : bool
                Enable or disable the channel

            Returns
            -------
            error_code : error_t
                0 for no error
                3 for invalid TC channel

            **/
            
            switch (_tc_pwm_channel) {
              case TC_PWM_CH_A: 
                if(enable)
                  TCCR4A |= (1<<COM4A1)|(1<<COM4A0)|(1<<PWM4A);
                else
                  // set bits 'COM4A1:0' to 0
                  TCCR4A &= ~((1<<COM4A1)|(1<<COM4A0)|(1<<PWM4A));
                break;
                
              case TC_PWM_CH_B: 
                if(enable)
                  TCCR4A |= (1<<COM4B1)|(1<<COM4B0)|(1<<PWM4B);
                else
                  // set bits 'COM4B1:0' to 0
                  TCCR4A &= ~((1<<COM4B1)|(1<<COM4B0)|(1<<PWM4B));
                break;

              case TC_PWM_CH_D: 
                if(enable)
                  TCCR4C |= (1<<COM4D1)|(1<<COM4D0)|(1<<PWM4D); 
                else
                  // set bits 'COM4D1:0' to 0
                  TCCR4C &= ~((1<<COM4D1)|(1<<COM4D0)|(1<<PWM4D)); 
                break;
              default: return InvalidChannel;
            }
            
            #if DEBUG
                if(enable)
                  Serial.print("Enable");
                else
                  Serial.print("Disable");
                Serial.print(" PWM output on TC channel ");
                Serial.println(_tc_pwm_channel);
            #endif
            
            return NoError;
        }
        
        error_t connect_pwm_output(byte io_pin) override {
            /**Connect internal Waveform Output 'OCnx' to physical 'io_pin'.
        
            Parameters
            ----------
            io_pin : byte
                Physical pin which will output PWM-generated waveform 'OCnx'
        
            Returns
            -------
            error_code : error_t
                0 for no error
                2 for invalid I/O pin
        
            **/
            
            #if DEBUG
              Serial.print("Connect PWM output pin ");
              Serial.println(io_pin);
            #endif
            
            switch(io_pin) {
              case 13: DDRC |= (1<<DDC7); break;
              case 10: DDRB |= (1<<DDB6); break;
              case  6: DDRD |= (1<<DDD7); break;
              default: return InvalidIOPin;
            }    
            return NoError;
        }     
        
        void enable_fast_pwm_mode() override {
            /**Enable Fast PWM Mode on TC4.
            TOP Counter stored in OCR4C, see Table 15-18.
            **/
        
            // set bits 'WGM11:0' to 0
            TCCR4D &= ~((1<<WGM41)|(1<<WGM40));
            TCCR4B |= (1<<PWM4X);
        }               

        error_t set_polarity_inverted(bool invert) override {
            /**Set the polarity of the PWM signal.
            Parameters
            ----------
            invert : bool
                Enable or disable inverted polarity
                
            Returns
            -------
            error_code : error_t
                0 for no error
                3 for invalid TC channel
            
            Notes
            -----
            Depending on the 'Compare Output Mode' settings COM4x1:0 (via
            Register TCCR4A and TCCR4C) the Output Compare pins OC4x will be
            set or cleared differently. Set COM4x1:0 to 0b10 to clear at match
            and set when TCNT4 = 0x000 (0b11: set at match, clear when 0x000).
            To avoid a voltage peak while changing the direction, due to
            change of _PWM_SIGN_PIN and PWM polarity are not in phase,
            set the counter TCNT4 to '_pwm_period-1' to reset it with the next
            tik.
            **/
            
            switch (_tc_pwm_channel) {
                case TC_PWM_CH_A:
                    if(!invert)
                        TCCR4A |= (1<<COM4A1)|(1<<COM4A0);
                    else {
                        TCCR4A |= (1<<COM4A1);
                        TCCR4A &= ~(1<<COM4A0);     // set bit 'COM4A0' to 0
                    } 
                    break;
                
                case TC_PWM_CH_B:
                    if(!invert)
                        TCCR4A |= (1<<COM4B1)|(1<<COM4B0);  
                    else {
                        TCCR4A |= (1<<COM4B1);
                        TCCR4A &= ~(1<<COM4B0);     // set bit 'COM4B0' to 0
                    }
                    break;
        
                case TC_PWM_CH_D:
                    if(!invert)
                        TCCR4C |= (1<<COM4D1)|(1<<COM4D0);
                    else {
                        TCCR4C |= (1<<COM4D1);
                        TCCR4C &= ~(1<<COM4D0);     // set bit 'COM4D0' to 0
                    }    
                    break;
                default: return InvalidChannel;
            }
            
            set_value(_pwm_period-1, _TCNTn);
            
            #if DEBUG
                if(invert)
                    Serial.print("Inverted ");
                else
                    Serial.print("Non-inverted ");
              Serial.print("polarity on TC channel ");
              Serial.println(_tc_pwm_channel);
            #endif
            
            return NoError;                                                               
        }
        
        void set_pwm_period(float frequency) {
            /**Set PWM period and PWM clock prescaler.
            
            Parameters
            ----------
            frequency : float
                PWM output frequency in Hz.
                
            Notes
            -----
            Notice that the resolution of the PWM scales with PWM frequency.
            Using prescalers for lower frequencies is recommended.
            For the 10-bit TC4 the clock prescaler value is set via pins
            CS43:CS40. The prescaler starts at 1 and with value 1 and doubles
            every step, i.e. 1->1, 2->2, 3->4, 4->8, 5->16, ...(see ATmega32U4
            datasheet table 15-14). Instead of dividing _pwm_period by 2
            a bitshift by 1 is more efficient.
            **/
            
            uint8_t clk_ps = 1;
            
            // see eq. (1) and (2) with no prescaler (clk_ps = 1)
            _pwm_period = (F_CPU_CLK / frequency) - 1;    // TOP = _pwm_period
            
            // If _pwm_period is larger than 10 bits, the prescaler is
            // increased and the _pwm_period divided by 2 (bitshift >> 1), to 
            // get the maximum possible resolution (see eq. (1), (2) and (3))
            while(_pwm_period & 0b1111110000000000) {
                clk_ps ++;
                _pwm_period = _pwm_period >> 1;
            }
            // To reset the prescaler (or counter TCNTn) set pins 
            // CS43:CS40 to zero (No clock source, Timer/Counter stopped)
            TCCR4B &= ~(0b1111<<CS40);        // set bits 'CS43:CS40' to 0
            TCCR4B |= clk_ps<<CS40;           // set prescaler

            set_value(_pwm_period, _TOP);     // TOP = _pwm_period
        }
        
        void set_duty_cycle_to_register(uint16_t value) override {
            set_value(value, _tc_reg_addr);
            
            #if DEBUG
              Serial.print("Set value ");
              Serial.print(value);
              Serial.print(" to register ");
              Serial.println(*_tc_reg_addr);
            #endif
        }
        
        // !!! USE AT YOUR OWN RISK !!!
        void enable_pwm_enhanced_mode() {
            /**Enable Enhanced PWM Mode on Timer/Counter4.

            Returns
            -------
            Nothing

            Notes
            -----
            In principle, this mode would allow to gain 1 more resolution
            bit for the PWM on Timer/Counter4 only. However, it seems like
            this mode has been NOT FUNCTIONING for Atmega32U4, at least up
            to revision D.
            **/
            // Enable Enhanced Compare/PWM mode (ENHC4 = 1)
            TCCR4E |= 1 << ENHC4;
        }
        
// -------------------------------- NOT USED --------------------------------
        
//         const uint16_t _TC_CLK_PRESCALER[15] = {1, 2, 4, 8, 16, 32, 64, 128,
//                                                 256, 512, 1024, 2048, 4096,
//                                                 8192, 16384};
        
//         uint16_t find_optimum_prescaler(float value) {
//             /**Find optimum clock prescaler.
//         
//             Parameters
//             ----------
//             value : float
//                 For calculating optimum clock prescaler
//         
//             Returns
//             -------
//             uint16_t : clk_ps
//                 Returns best prescaler
//         
//             Notes
//             -----   
//             Calculate the TC4 clock prescaler: Use the smallest possible
//             prescaler which will result in the highest achievable PWM
//             resolution, i.e. make OCR4C as large as possible.
//             For a given CPU and PWM frequency, choose a clock prescaler that
//             keeps the PWM resolution as close as possible to the maximum
//             of 10 bits.
//             For example, a 16 MHz CPU clock and 20 kHz PWM yields 800 steps
//             in the interval [0,799]. This corresponds to a resolution of
//             9.64 bits.
//             **/
//             
//             float pwm_resolution;
//             float f_clk_tc;
//             uint16_t clk_ps;
// 
// 
//             // Iterate through the available TC1 prescaler values:
//             for(uint8_t i = 0; i < 15; i++) {
//                       
//                 clk_ps = _TC_CLK_PRESCALER[i];                
//                 f_clk_tc = F_CPU_CLK / clk_ps;              // eq. (2)
//                 pwm_resolution = log2(f_clk_tc / value);    // eq. (3)
//                 
//                 // Jump out of the loop, if best prescaler value is found.
//                 if(pwm_resolution <= 10.0) break;
//             }
//             
//             #if DEBUG
//               Serial.print("Found ");
//               Serial.print(clk_ps);
//               Serial.print(" as optimum prescaler for pwm frequency ");
//               Serial.println(value);
//             #endif
//             
//             return clk_ps;
//         }
// 
//         error_t set_pwm_clock_prescaler(uint16_t tc_clock_prescaler) {
//             /**Set PWM clock prescaler.
//         
//             Parameters
//             ----------
//             tc_clock_prescaler : uint16_t
//                 Sets the requested TC clock prescaler
//         
//             Returns
//             -------
//             error_code : error_t
//                 0 for no error
//                 4 for invalid clock prescaler
//             
//             **/
//             
//             byte ps_select;
//         
//             // Convert the TC4 clock prescaler value to the respective 
//             // prescaler select entry, see ATmega32U4 datasheet table 15-14
//             switch(tc_clock_prescaler) {
//               case     1: ps_select = 1; break;
//               case     2: ps_select = 2; break;
//               case     4: ps_select = 3; break;
//               case     8: ps_select = 4; break;
//               case    16: ps_select = 5; break;
//               case    32: ps_select = 6; break;
//               case    64: ps_select = 7; break;
//               case   128: ps_select = 8; break;
//               case   256: ps_select = 9; break;
//               case   512: ps_select = 10; break;
//               case  1024: ps_select = 11; break;
//               case  2048: ps_select = 12; break;
//               case  4096: ps_select = 13; break;
//               case  8192: ps_select = 14; break;
//               case 16384: ps_select = 15; break;
//               default: return InvalidPrescaler;
//             }
//             
//             // To reset the prescaler (or counter TCNTn) set pins 
//             // CS43:CS40 to zero (No clock source, Timer/Counter stopped)
//             TCCR4B &= ~(0b1111<<CS40);        // set bits 'CS43:CS40' to 0
//             TCCR4B |= ps_select<<CS40;        // set prescaler
//             
//             #if DEBUG
//               Serial.print("Write ");
//               Serial.print(ps_select,BIN);
//               Serial.print(" to register TCCR4B for prescaler ");
//               Serial.println(tc_clock_prescaler);
//             #endif
//             
//             return NoError;          
//         }

// --------------------------------------------------------------------------
                
};

// -------------------------------- NOT USED --------------------------------

// void tc3_sfr_reset() {
// Some special function registers belonging to
// Timer/Counter3 are not properly zero'ed on POR
//     TCCR3A = 0;
//     TCCR3B = 0;
//     TCCR3C = 0;
//     OCR3AH = 0;
//     OCR3AL = 0;
//     OCR3BH = 0;
//     OCR3BL = 0;
//     OCR3CH = 0;
//     OCR3CL = 0;
//     TIMSK3 = 0;
//     TIFR3 = 0;
// }

// --------------------------------------------------------------------------

#endif  // _VIPERLEED_B_FIELD_PWM