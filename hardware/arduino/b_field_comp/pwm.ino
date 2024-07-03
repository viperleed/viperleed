/*
ViPErLEED - Driver for using fast-PWM mode via Timer/Counter4.
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/

#include <Arduino.h>        // For interrupts()/noInterrupts()
#include "pwm.h"


/**
Quick Start Guide to Timer/Counter4
===================================
This section is intended as a concise introduction to the PWM waveform
generation using the Timer/Counter4, subsequently abbreviated as TC4.

TC4 can be used in 8-bit or 10-bit operation, depending on whether the user 
also programs the two high bits inside shared register TC4H (cf. datasheet 
sec. 15.2.2). For example, in order to perform the 10-bit write 0x31F to 
register TC4H:OCR4C, first program TC4H with 0x3, then write 0x1F to OCR4C.

TC4 contains three independant Output Compare Units A/B/D which constantly
compare their register values OCR4A/B/D with monotonically increasing or
decreasing Counter register TCNT4. If a match occurs, Output Compare pins
OC4A/B/D will be set or cleared depending on the Compare Output Mode settings 
COM4x1:0, where x = Output Compare Unit A/B/D. COM4x1:0 can thus be used to 
invert the  PWM signal, see 'set_pwm_polarity'. Also, the COM4x1:0 bits control 
the connection to the inverted/non-inverted Output Compare pins /OC4x vs OC4x.

While registers OCR4A/B/D are used to control the PWM duty cycle, OCR4C is used
to set the PWM base frequency, and implicitly also the PWM resolution. The 
relevant formulae are also given in the datasheet sec. 15.8. For TC4 operated
in Fast PWM Mode the PWM frequency, TC4 clock frequency and PWM resolution are
determined as follows:

(1)   f_PWM = f_clk_T4 / (OCR4C + 1)

(2)   f_clk_T4 = F_CPU_CLK / TC4_prescaler

(3)   res_PWM = log2(OCR4C + 1) = log2(f_clk_T4 / f_PWM)

Every time TCNT4 has a match with OCR4C, TCNT4 is reset to zero and begins to
count up or down again until another match occurs. The TC4 frequency 'f_clk_T4'
is derived from the CPU clock 'F_CPU_CLK' divided by the TC4 prescaler settings
CS43:40 contained in register TCCR4B. Equation (3) shows that the user may want 
to stay as close to the maximum value of OCR4C as possible, in order to achieve 
the highest possible PWM resolution. Therefore, 'set_pwm_frequency' will first
determine the TC4 prescaler value which will allow for the highest possible
10-bit value inside TC4H:OCR4C which in turn will provide the best possible PWM
resolution. If a value smaller than three is written to OCR4C, the hardware will
replace that value with three in order to prevent too high a PWM frequency.

The ATmega32U4 includes the following Timer/Counter modules:
- 8-bit Timer/Counter0
- 16-bit Timers/Counters (TC1 and TC3)
- 10-bit High Speed Timer/Counter4

In contrast to the other TC modules, TC4 can be operated with up to 64 MHz if
the supply voltage is higher than 4V. The firmware does not make use of this 
feature.

Note: The 11-bit Enhanced PWM Mode is not working up to some chip revision.
**/

error_t set_pwm_frequency(double f_pwm) {
    /**Set PWM frequency, using 10-bit Timer/Counter4.

    Parameters
    ----------
    f_pwm : double
        PWM output frequency in Hz.
        Notice that the resolution of the PWM scales with PWM frequency.
        The PWM frequency is shared across all three available Timer/Counter4
        channels A, B and D.

    Returns
    -------
    error_code : byte
        0 for no error
        1 for requested pwm frequency out-of-range

    **/

    double f_clk_t4;
    uint16_t clk_ps;
  
    if (f_pwm < F_PWM_MIN || f_pwm > F_PWM_MAX)
        return OutOfRange;

    enable_fast_pwm_mode();  
 
    clk_ps = find_optimum_prescaler(f_pwm);
    set_pwm_clock_prescaler(clk_ps);

    // Do not use "enhanced mode". Functional only for certain chip revisions.
    // enable_pwm_enhanced_mode();

    f_clk_t4 = F_CPU_CLK / clk_ps;

    // The PWM period (= 1 / f_pwm) is defined in terms of TC4 clock ticks:
    pwm_period = f_clk_t4 / f_pwm - 1;
    set_ten_bit_value(pwm_period, &OCR4C);

    connect_pwm_output(COIL_1_PWM);
    connect_pwm_output(COIL_2_PWM);

    return NoError;
}


error_t set_signed_pwm_value(double value, byte sign_select_pin, byte *tc4_reg_addr) {
    /**Set PWM duty cycle (i.e., average voltage), including sign output.

    Parameters
    ----------
    value : double
        PWM duty cycle to be set. Should be between -1.0 and +1.0.
    sign_select_pin : byte
        The Arduino pin that takes care of the sign of this PWM signal.
    tc4_reg_addr : byte*
        Address of the specific Timer/Counter4 register.

    Returns
    -------
    error_code : byte
        0 for no error
        1 for 'value' (aka coil current) out-of-range
    **/

    if (value < -1.0 || value > 1.0) 
        return OutOfRange;
 
    value = set_current_sign(value, sign_select_pin);
     
    // Notice that the 'coil_current', i.e. the time-averaged value of the signal
    // from the PWM, is exactly the same as the duty cycle of the PWM itself.
    // First, set or clear an I/O pin depending on the current direction;
    // Then set the PWM duty cycle for 'COIL_1' or 'COIL_2' at 'tc4_reg_addr':
    // The duty cycle is the desired 'on' time in percent, i.e. 'coil_current'
    // times the PWM period measured in TC4 clock ticks.
    // Also note it is sufficient to simply pass 'abs(value)' to program the 
    // duty cycle because we set the current direction with 'set_pwm_polarity'.
    set_ten_bit_value((uint16_t)round(value * (pwm_period)), tc4_reg_addr);

    return 0;
}


double set_current_sign(double value, byte sign_select_pin) {
    /**Set digital I/O pin according to current direction in a coil.

    Parameters
    ----------
    value : double
        Coil current direction: Either positive or negative
    sign_select_pin : byte
        Which I/O pin to use as a sign indicator
    **/
    if (value < 0) {
        digitalWrite(sign_select_pin, HIGH);
        return value + 1;
    } 
    digitalWrite(sign_select_pin, LOW);
    return value;
}


void set_ten_bit_value(uint16_t ten_bit_value, volatile uint8_t *reg) {
    /**Write 10-bit value to Timer/Counter4 register.

    Parameters
    ----------
    ten_bit_value : uint16_t
        10-bit value to be written
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

    noInterrupts();  // Ensure nothing bothers setting two registers
    TC4H = ten_bit_value >> 8;
    *reg = ten_bit_value & 255;
    interrupts(); 
}


error_t connect_pwm_output(byte io_pin) {
    /**Connect internal Waveform Output 'OCW4x' to 'io_pin'.

    Parameters
    ----------
    io_pin : byte
        Physical pin which will output PWM-generated waveform 'OCW4x'

    Returns
    -------
    error_code : error_t
        0 for no error
        2 for invalid I/O pin

    Notes
    -----
    Pins OC4A/B/D: cleared on compare match (TCNT4 == OCR4A/B/D), set when
    TCNT4 = 0x000; Connect the Waveform Outputs OCW4A/B/D to Output Compare
    pins OC4A/B/D. For more info see Atmega32U4 datasheet, section 15.12.1
    **/
    
    // Route the Timer/Counter4 Waveform Outputs OCW4A/B/D to the non-inverting
    // Output Compare pins OC4A/B/D: 
    switch(io_pin) {
      case  6: TCCR4C |= (1 << COM4D1) | (0 << COM4D0); break;
      case 10: TCCR4A |= (1 << COM4B1) | (0 << COM4B0); break;
      case 13: TCCR4A |= (1 << COM4A1) | (0 << COM4A0); break;
      default: return InvalidIOPin;
    }    
    return NoError;
}


  // Calculate the TC4 clock prescaler: Use the smallest possible prescaler
  // which will result in the highest achievable PWM resolution, i.e. make
  // OCR4C as large as possible. Note: Some TC4 registers including OCR4C are
  // extended to 10 bits by adding two bits in TC4H (datasheet sec. 15.2.2).
  // Relevant formulae: 
  // f_OC4x_pwm = f_clk_T4 / (N + 1), where N ... value of OCR4C
  // res_pwm = log2(OCR4C + 1), in bits. See datasheet sec. 15.8.
  uint16_t find_optimum_prescaler(double f_pwm) {
    double pwm_resolution;
    double f_clk_t4;
    uint16_t clk_ps;

    // For a given CPU and PWM frequency, choose a clock prescaler that keeps
    // the PWM resolution as close as possible to the maximum of 10 bits.
    // Iterate through the available TC4 prescaler values:
    for(uint8_t i = 0; i < 15; i++)
    {
      f_clk_t4 = F_CPU_CLK / TC4_CLK_PRESCALER[i];
      pwm_resolution = log2(f_clk_t4 / f_pwm);
      clk_ps = TC4_CLK_PRESCALER[i];

      // Jump out of the loop if we have found the best prescaler value
      if(pwm_resolution <= 10.0) break;
    }
    return clk_ps;
  }


error_t set_pwm_clock_prescaler(uint16_t tc4_clock_prescaler) {
    /**Set PWM clock prescaler.

    Parameters
    ----------
    tc4_clock_prescaler : uint16_t
        Sets the requested TC4 clock prescaler

    Returns
    -------
    error_code : error_t
        0 for no error
        4 for invalid clock prescaler

    Notes
    -----
    TC4 clock frequency = CPU clock frequency / TC4 clock prescaler.
    Notice that we should divide the CPU clock by the smallest amount 
    possible: the faster the TC4 counter, the better the PWM resolution.
    For example, a 16 MHz CPU clock and 20 kHz PWM yields 800 steps in
    the interval [0,799]. This corresponds to a resolution of 9.64 bits.
    **/
    byte ps_select;

 
    // Convert the TC4 clock prescaler value to the respective 
    // prescaler select entry, see ATmega32U4 datasheet table 15-14
    switch(tc4_clock_prescaler) {
      case 1: ps_select = 1; break;
      case 2: ps_select = 2; break;
      case 4: ps_select = 3; break;
      case 8: ps_select = 4; break;
      case 16: ps_select = 5; break;
      case 32: ps_select = 6; break;
      case 64: ps_select = 7; break;
      case 128: ps_select = 8; break;
      case 256: ps_select = 9; break;
      case 512: ps_select = 10; break;
      case 1024: ps_select = 11; break;
      case 2048: ps_select = 12; break;
      case 4096: ps_select = 13; break;
      case 8192: ps_select = 14; break;
      case 16384: ps_select = 15; break;  
      default: return InvalidPrescaler;
    }
    TCCR4B |= ps_select;
    return NoError;
}


error_t enable_pwm_channel(TC4_PWM_CHANNEL channel, bool enable) {
    /**Enable/disable PWM output on TC4 channels A, B or D.

    Parameters
    ----------
    channel : TC4_PWM_CHANNEL
        Selects the TC4 channel to enable/disable
    enable : bool
        Enable or disable the channel

    Returns
    -------
    error_code : error_t
        0 for no error
        3 for invalid TC4 channel

    Notes
    -----
    The register values of OCR4A/B/D will determine
    the duty cycle of the PWM: When the TC4 counter 
    register TCNT4 reaches the values in these registers,
    pins OC4A/B/D will toggle.
    More info: see Atmega32U4 datasheet, section 15.8.2
    **/
    switch (channel) {
      case TC4_PWM_CH_A: 
        if(enable)
          TCCR4A |= 1 << PWM4A; 
        else
          TCCR4A &= ~(1 << PWM4A);
        break;

      case TC4_PWM_CH_B: 
        if(enable)
          TCCR4A |= 1 << PWM4B;
        else
          TCCR4A &= ~(1 << PWM4B);    
        break;

      case TC4_PWM_CH_D: 
        if(enable)
          TCCR4C |= 1 << PWM4D; 
        else
          TCCR4C &= ~(1 << PWM4D); 
        break;
      default: return InvalidChannel;
    }
    return NoError;
}


void enable_fast_pwm_mode() {
    /**Enable Fast PWM Mode on TC4.

    Returns
    -------
    Nothing
    **/
    // Keep PWM output non-inverted (Clear 'PWM4X')
    TCCR4B &= ~(1 << PWM4X);

    // Enable Fast PWM Mode (WGM41:40 = 0B00)
    TCCR4D &= ~((1 << WGM41) | (1 << WGM40));  // Clear <WGM41:WGM40> only
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


uint8_t pin_to_tc4_reg_addr(uint8_t pwm_pin) {
    switch(pwm_pin) {
      case  6: return _SFR_ADDR(OCR4D);
      case 10: return _SFR_ADDR(OCR4B);
      case 13: return _SFR_ADDR(OCR4A);
      default: return 0;                                                       // Change to more descriptive error
    }
}


// If pin definitions on the ATmega32U4 should change, 
// the assignments below would have to change accordingly.
TC4_PWM_CHANNEL pin_to_tc4_channel(uint8_t pwm_pin) {
    switch(pwm_pin) {
      case  6: return TC4_PWM_CH_D;
      case 10: return TC4_PWM_CH_B;
      case 13: return TC4_PWM_CH_A;
      default: return 0;                                                       // Change to more descriptive error
    }
}


// Some special function registers belonging to
// Timer/Counter4 are not properly zero'ed on POR
void tc4_sfr_reset() {
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


double log2(double val)
{
   return log(val) / log(2);
}
