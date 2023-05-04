/**                                                                             // TODO: add description of the contents of the file + autors/date...
**/

#include <Arduino.h>  // for interrupts()/noInterrupts()


#define COIL_1  6   // Also A7;  PD7 on Atmega32U4
#define COIL_2 10   // Also A10; PB6 on Atmega32U4

// Current direction: positive or negative?
#define COIL_1_SIGN 18  // PF7 on ATmega32U4; also used for INA on shunt
#define COIL_2_SIGN 19  // PF6 on ATmega32U4; also used for INA on shunt
#define POSITIVE_CURRENT  1
#define NEGATIVE_CURRENT -1

#define F_CLK_T4     16000      // Timer/Counter4 clock = 16 MHz
#define PWM_MIN_FREQ 15.625     // Value comes from maximum value for
                                // register OCR4C and from the choice of a
                                // value of 1 in set_pwm_clock_prescaler


// pwm_clock_divider: how many intervals the counters count before
// rolling over. This essentially sets the frequency of the PWM.
// 'uint16_t' is sufficient because TC4H:OCR4C is 10 bits wide.
uint16_t pwm_clock_divider;  // Use set_pwm_frequency for setting!


void setup()
{
    set_pwm_frequency(20);  // 20 kHz

    set_coil_current(0.625, COIL_1);
    set_coil_current(0.25,  COIL_2);

    pinMode(COIL_1, OUTPUT);             // Define PD7 (OC4D) as output
    pinMode(COIL_2, OUTPUT);             // Define PB6 (OC4B) as output

    pinMode(COIL_1_SIGN, OUTPUT);     // Define PF7 as output
    pinMode(COIL_2_SIGN, OUTPUT);     // Define PF6 as output
}


void loop()
{
    delay(1000);
}


byte set_pwm_frequency(double freq){
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

    // How many times freq fits in F_CLK_T4. Notice the -1: the PWM
    // counter will reset when it counts 'pwm_clock_divider+1' intervals
    pwm_clock_divider = F_CLK_T4 / freq - 1;

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
    // E.g. OCR4A = Output Compare register on Timer/Counter 4, channel A

    // Set PWM frequency to freq (TC4H:OCR4C = 799 gives freq exactly;
    // datasheet formula [Section 15.8.2] off by one for fast PWM mode)
    set_ten_bit_value(pwm_clock_divider, &OCR4C);
}


byte set_coil_current(double coil_current, uint8_t coil){
    /**Set fraction of maximum current in a coil.

    Parameters
    ----------
    coil_current : double
        Fraction of maximum current. Should be between zero and one.
    coil : {COIL_1, COIL_2}
        Which coil's current should be set.

    Returns
    -------
    error_code : byte
        0 for no error
        2 for coil_current out-of-range
        3 for invalid coil
    **/
    uint8_t *_reg_addr;
    byte sign_select_pin;

    if (coil_current < -1 || coil_current > 1) return 2;                        // TODO: could make these return values into error codes, similar to the driver codes, or use a bunch of defines
    switch(coil)
    {
        case COIL_1:
            _reg_addr = &OCR4D;
            sign_select_pin = COIL_1_SIGN;
            break;
        case COIL_2:
            _reg_addr = &OCR4B;
            sign_select_pin = COIL_2_SIGN;
            break;
        default:
            return 3;
    }

    // Notice that the coil_current, i.e., the time-averaged
    // value of the signal from the PWM, is exactly the same
    // as the duty cycle of the PWM itself.
    // Set PWM duty cycle for channel at `_reg_addr`:
    //    `duty_cycle` == `coil_current` = TC4H:`_register` / TC4H:OCR4C
    set_current_sign(coil_current, sign_select_pin);
    set_ten_bit_value(pwm_clock_divider * coil_current, _reg_addr);
    return 0;
}


void set_current_sign(byte sign, byte sign_select_pin){
    if (sign < 0){
        set_pwm_polarity(NEGATIVE_CURRENT);
        digitalWrite(sign_select_pin, HIGH);                                    // TODO: check if this is the right way (test on coil)
    }
    else {
        set_pwm_polarity(POSITIVE_CURRENT);
        digitalWrite(sign_select_pin, LOW);
    }
}


void set_ten_bit_value(uint16_t ten_bits_value, uint8_t *REGISTER){
    // Registers are 8 bits. Some can also accept 10-bit values.
    // To do this, the two highest bits are to be written in the
    // *shared* register TC4H right before the remaining 8 bits
    // are written in the desired register. More info: Atmega32U4
    // datasheet, section 15.11.
    noInterrupts();  // Ensure nothing bothers setting two registers
    TC4H = ten_bits_value >> 8;
    *REGISTER = ten_bits_value & 255;
    interrupts();
}


void set_pwm_polarity(byte polarity){
    // Pins OC4B, OC4BD: cleared on compare match (TCNT = OCR4B/D),
    // set when TCNT = 0x000; Enable PWM output channels B and D
    // This part essentially selects whether we output "high" or
    // "low" when the counter reaches the threshold. This can be
    // used for flipping the signal.
    // For more info, see Atmega32U4 datasheet, section 15.12.1
    TCCR4C &= ~((1 << COM4D1) | (1 << COM4D0));   // Clear <COM4D1:COM4D0>
    TCCR4A &= ~((1 << COM4B1) | (1 << COM4B0));   // Clear <COM4B1:COM4B0>

    // Set duty cycle polarity: for positive currents, we reset to
    // zero when the count goes above the threshold; for negative
    // currents, we set to one when the count goes above threshold
    if (polarity == POSITIVE_CURRENT)
        polarity = 0;
    else
        polarity = 1;

    // For both channels B and D, activate only the non-inverting
    // output; set where switch occurs depending on polarity
    TCCR4A |= (1 << COM4B1) | (polarity << COM4B0);
    TCCR4C |= (1 << COM4D1) | (polarity << COM4D0);
}


void set_pwm_clock_prescaler(){
    // Set Timer/Counter4 prescaler to 1;
    // Timer/Counter4 clock frequency = system clock frequency / 1
    // (f_clk_T4 = 16 MHz)
    // Notice that we should not divide it further: the faster the counter,
    // the better our resolution. With 16 MHz input and 20 kHz output one
    // gets ~ 16 MHz / 20 kHz = 800 steps of resolution.

    // NOTE: Should one decide to use a different clock divider, the value
    //       of the PWM_MIN_FREQ should be changed accordingly!
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
    // Make OCR4B/D the registers whose values will be
    // used to determine the duty cycle of the PWM. I.e.,
    // when the corresponding counter TCNT4 reaches the
    // values in these registers pins OC4B/D will toggle
    // More info: see Atmega32U4 datasheet, section 15.8.2
    TCCR4A |= (1 << PWM4B);
    TCCR4C |= (1 << PWM4D);
}


void set_fast_pwm_mode(){
    // Keep PWM output non-inverted (Clear 'PWM4X')
    TCCR4B &= ~(1 << PWM4X);

    // Enable Fast PWM Mode (WGM41:40 = 0B00)
    TCCR4D &= ~((1 << WGM41) | (1 << WGM40));  // Clear <WGM41:WGM40> only
}


// !!! USE AT YOUR OWN RISK !!!
void use_pwm_enhanced_mode(){
    // In principle, this mode would allow to gain 1 more resolution
    // bit for the PWM on Timer/Counter4 only. However, it seems like
    // this mode has been NOT FUNCTIONING for Atmega32U4, at least up
    // to revision D.
    // Enable Enhanced Compare/PWM mode (ENHC4 = 1)
    TCCR4E |= 1 << ENHC4;
}
