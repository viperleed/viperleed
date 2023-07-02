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
    freq : double
        PWM output frequency in kilohertz. Minimum: 15.625 Hz.
        Notice that the resolution of the PWM (i.e., of the
        average current) scales with freq. There can only be
        F_CLK_T4 / freq - 1 individual steps.

    Returns
    -------
    error_code : byte
        0 for no-error
        1 freq out-of-range

    **/
    if (freq < PWM_MIN_FREQ || freq > F_CLK_T4) return 1;


    set_pwm_polarity(POSITIVE_CURRENT);
    set_pwm_threshold_channels();
    set_fast_pwm_mode();
    set_pwm_clock_prescaler();

    // Don't use "enhanced mode". Seems functional
    // only for certain chip revisions.
    // use_pwm_enhanced_mode();

    // Atmega32U4 register names 'OCRnx' contain device number 'n'
    // (where 'n' is Timer/Counter n) and Output Compare unit 'x'
    // (where 'x' is A/B/C)
    // E.g. OCR4A = Output Compare register on Timer/Counter4, channel A

    // Set PWM frequency to freq (TC4H:OCR4C = 799 gives freq exactly;
    // datasheet formula [Section 15.8.2] off by one for fast PWM mode)
    return 0;
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


void set_ten_bit_value(uint16_t ten_bits_value, volatile uint8_t *REGISTER){
    /**Write 10-bit value to Timer/Counter4 register.

    Parameters
    ----------
    ten_bits_value : uint16_t
    REGISTER : uint8_t *
        10-bit value to be written
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
    TC4H = ten_bits_value >> 8;
    *REGISTER = ten_bits_value & 255;
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
    TCCR4A |= (1 << COM4B1) | (polarity << COM4B0);
    TCCR4C |= (1 << COM4D1) | (polarity << COM4D0);
}


void set_pwm_clock_prescaler(){
    /**Set PWM clock prescaler.

    Returns
    -------
    Nothing

    Notes
    -----
    - Set Timer/Counter4 prescaler to 1;
    - TC4 clock frequency = CPU clock frequency / TC4 clock prescaler
      Notice that we should divide the CPU clock by the smallest amount 
      possible: the faster the TC4 counter, the better the PWM resolution.
      For example, with a 16 MHz CPU clock and 20 kHz output we get 16 MHz
      divided by 20 kHz = 800 steps of resolution (9.64 bits)

    **/
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
    TCCR4B |= (0 << CS43) | (0 << CS42) | (0 << CS41) | (1 << CS40);
}


void set_pwm_threshold_channels(){
    /**Enable PWM output on Timer/Counter4 channels B and D.

    Returns
    -------
    Nothing

    Notes
    -----
    OCR4A/B/D are the registers whose values will be
    used to determine the duty cycle of the PWM. I.e.,
    when the TC4 counter register TCNT4 reaches the
    values in these registers pins OC4A/B/D will toggle
    More info: see Atmega32U4 datasheet, section 15.8.2
    **/
    TCCR4A |= (1 << PWM4B);
    TCCR4C |= (1 << PWM4D);
}


void set_fast_pwm_mode(){
    /**Enable Fast PWM Mode on Timer/Counter4.

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
