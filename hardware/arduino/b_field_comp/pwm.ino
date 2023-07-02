/*
ViPErLEED - Driver for using fast-PWM mode via Timer/Counter4.
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/

#include <Arduino.h>        // For interrupts()/noInterrupts()
#include "pwm.h"


byte set_pwm_frequency(double f_pwm) {
    /**Set PWM frequency, using the 10-bit Timer/Counter4.

    Parameters
    ----------
    f_pwm : double
        PWM output frequency in Hz.
        Notice that the resolution of the PWM (i.e., of the average current)
        scales with PWM frequency. The PWM frequency is shared across all 
        three available Timer/Counter4 channels A, B and D.

    Returns
    -------
    error_code : byte
        0 for no-error
        1 freq out-of-range

    **/
  double pwm_resolution;
  double f_clk_t4;
  uint16_t clk_ps;
  
  if (f_pwm < F_PWM_MIN || f_pwm > F_PWM_MAX) return 1;

  // Calculate the TC4 clock prescaler: Use the smallest possible prescaler
  // which will result in the highest achievable PWM resolution, i.e. make
  // OCR4C as large as possible. Note: Some TC4 registers such as OCR4C are
  // extended to 10 bits by adding two bits in TC4H (datasheet sec. 15.2.2).
  // Relevant formulae: 
  // f_OC4nx_pwm = f_clk_T4 / (N + 1), where N ... value of OCR4C
  // res_pwm = log2(OCR4C + 1), in bits. See datasheet sec. 15.8.
  
  // For a given CPU and PWM frequency, choose a clock prescaler that keeps
  // the PWM resolution as close as possible to the maximum of 10 bits.
  // Iterate through the available TC4 prescaler values:
  for(uint8_t i = 0; i < 16; i++)
  {
    f_clk_t4 = F_CPU_CLK / TC4_CLK_PRESCALER[i];
    pwm_resolution = log2(f_clk_t4 / f_pwm);
    clk_ps = TC4_CLK_PRESCALER[i];

    //Serial.print(clk_ps,DEC); Serial.print(": "); Serial.println(pwm_resolution,3);

    // Jump out of the loop if we have found the best prescaler value
    if(pwm_resolution <= 10.0) break;
  }
  set_pwm_clock_prescaler(clk_ps);

  // Don't use "enhanced mode". Functional only for certain chip revisions.
  // use_pwm_enhanced_mode();

  // The PWM period (= 1 / f_pwm) is defined in terms of TC4 clock ticks:
  pwm_period = f_clk_t4 / f_pwm - 1;  
  set_ten_bit_value(pwm_period, &OCR4C);
}


byte set_signed_pwm_value(double value, byte sign_select_pin, byte *tc4_reg_addr) {
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
        2 for coil_current out-of-range
    **/
    if (value < -1.0 || value > 1.0) return 2;                                      // TODO: could make these return values into error codes, similar to the driver codes, or use a bunch of defines

  // Notice that the 'coil_current', i.e. the time-averaged value of the signal
  // from the PWM, is exactly the same as the duty cycle of the PWM itself.
  // First, set or clear an I/O pin depending on the current direction;
  // Then set the PWM duty cycle for 'COIL_1' or 'COIL_2' at 'tc4_reg_addr':
  // The duty cycle is the desired ON time in percent, i.e. 'coil_current'
  // times the PWM period measured in TC4 clock ticks.
  set_current_sign(value, sign_select_pin);
  set_ten_bit_value(value * pwm_period, tc4_reg_addr);
  return 0;
}



void set_current_sign(double value, byte sign_select_pin) {
    /**Set digital I/O pin according to current direction in a coil.

    Parameters
    ----------
    value : double
        Coil current direction: Either positive or negative
    sign_select_pin : byte
        Which I/O pin to use as a sign indicator
    **/
    if (value < 0){
        set_pwm_polarity(NEGATIVE_CURRENT);
        digitalWrite(sign_select_pin, HIGH);                                    // TODO: check if this is the right way (test on coil)
    }
    else {
        set_pwm_polarity(POSITIVE_CURRENT);
        digitalWrite(sign_select_pin, LOW);
    }
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


void set_pwm_polarity(byte polarity) {
    /**Set PWM polarity.

    Parameters
    ----------
    polarity : byte
        Sets the requested counter polarity

    Notes
    -----
    Pins OC4A/B/D: cleared on compare match (TCNT4 == OCR4A/B/D), set when
    TCNT4 = 0x000; Connect the Waveform Outputs OCW4A/B/D to Output Compare
    pins OC4A/B/D.
    This function essentially selects whether we output "high" or "low" when
    TCNT4 reaches the threshold. This can be used for inverting the signal.
    For more info see Atmega32U4 datasheet, section 15.12.1
    **/
    TCCR4A &= ~((1 << COM4A1) | (1 << COM4A0));     // Clear <COM4A1:COM4A0>
    TCCR4A &= ~((1 << COM4B1) | (1 << COM4B0));     // Clear <COM4B1:COM4B0>
    TCCR4C &= ~((1 << COM4D1) | (1 << COM4D0));     // Clear <COM4D1:COM4D0>

    // Set duty cycle polarity: for positive currents, we reset to
    // zero when the count goes above the threshold; for negative
    // currents, we set to one when the count goes above threshold
    if (polarity == POSITIVE_CURRENT)
        polarity = 0;
    else
        polarity = 1;

    // For channels A, B and D, route the Waveform Outputs OCW4A/B/D to the
    // non-inverting Output Compare pins OC4A/B/D; 
    // Set where pin toggle occurs depending on polarity:
    TCCR4A |= (1 << COM4A1) | (polarity << COM4A0);
    TCCR4A |= (1 << COM4B1) | (polarity << COM4B0);
    TCCR4C |= (1 << COM4D1) | (polarity << COM4D0);
}


void set_pwm_clock_prescaler(uint16_t tc4_clock_prescaler) {
    /**Set PWM clock prescaler.

    Parameters
    ----------
    tc4_clock_prescaler : uint16_t
        Sets the requested TC4 clock prescaler

    Notes
    -----
    - Set Timer/Counter4 prescaler to 1;
    - TC4 clock frequency = CPU clock frequency / TC4 clock prescaler
      Notice that we should divide the CPU clock by the smallest amount 
      possible: the faster the TC4 counter, the better the PWM resolution.
      For example, with a 16 MHz CPU clock and 20 kHz output we get 16 MHz
      divided by 20 kHz = 800 steps of resolution (9.64 bits)

    **/
    uint8_t ps_select;
    // NOTE: Power-on-reset should zero the entire register, but it seems
    // some chip revisions don't do that properly.
    // Clear entire register except MSbit 'PWM4X'
    TCCR4B &= ~((1 << PSR4)
                | (1 << DTPS41)
                | (1 << DTPS40)
                | (1 << CS43)
                | (1 << CS42)
                | (1 << CS41)
                | (1 << CS40));

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
    }
    TCCR4B |= ps_select;
}


void enable_pwm_channel(TC4_PWM_CHANNEL channel, bool enable) {
    /**Enable/disable PWM output on TC4 channels A, B or D.

    Parameters
    ----------
    channel : TC4_PWM_CHANNEL
        Selects the TC4 channel to enable/disable
    enable : bool
        Enable or disable the channel

    Notes
    -----
    OCR4A/B/D are the registers whose values will be
    used to determine the duty cycle of the PWM. I.e.,
    when the TC4 counter register TCNT4 reaches the
    values in these registers pins OC4A/B/D will toggle
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
    }
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
void use_pwm_enhanced_mode() {
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
// TODO: Move to somewhere else?
double log2(double val)
{
   return log(val) / log(2);
}
