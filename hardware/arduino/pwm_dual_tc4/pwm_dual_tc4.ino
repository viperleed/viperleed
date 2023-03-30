

#define COIL_1 PD7   // Also A7,  and D6
#define COIL_2 PB6   // Also A10, and D10

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

  pinMode(COIL_1, OUTPUT);             // Define PC7 (OC4A) as output
  pinMode(COIL_2, OUTPUT);             // Define PC7 (OC4A) as output
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

  set_pwm_polarity();
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
  set_ten_bit_value(pwm_clock_divider, OCR4C);
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
  uint8_t _register;

  if (coil_current < 0 || coil_current > 1) return 2;
  switch(coil)
  {
    case COIL_1:
      _register = OCR4D;
      break;
    case COIL_2:
      _register = OCR4B;
      break;
    default:
      return 3;
  }

  // Notice that the coil_current, i.e., the time-averaged
  // value of the signal from the PWM, is exactly the same
  // as the duty cycle of the PWM itself.
  // Set PWM duty cycle for channel `_register`:
  //    `duty_cycle` == `coil_current` = TC4H:`_register` / TC4H:OCR4C
  set_ten_bit_value(pwm_clock_divider  * coil_current, _register);
  return 0;
}


void set_ten_bit_value(uint16_t ten_bits_value, uint8_t REGISTER){
  // Registers are 8 bits. Some can also accept 10-bit values.
  // To do this, the two highest bits are to be written in the
  // *shared* register TC4H right before the remaining 8 bits
  // are written in the desired register. More info: Atmega32U4
  // datasheet, section 15.11.
  TC4H = ten_bits_value >> 8;
  REGISTER = ten_bits_value & 255;
}


void set_pwm_polarity(){
  // Pins OC4A, OC4B: cleared on compare match (TCNT = OCR4A/B),
  // set when TCNT = 0x000; Enable PWM output channels A and B
  // This part essentially selects whether we output "high" or
  // "low" when the counter reaches the threshold. This can be
  // used for flipping the signal.
  // For more info, see Atmega32U4 datasheet, section 15.12.1
  TCCR4A |= (1 << COM4A1)   | (0 << COM4A0)  // count up, non-inverting, OC4A
            | (1 << COM4B1) | (0 << COM4B0); // count up, non-inverting, OC4B
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
  TCCR4B |= (0 << CS43) | (0 << CS42) | (0 << CS41) | (1 << CS40);
}


void set_pwm_threshold_channels(){
  // Make OCR4A/B the registers whose values will be
  // used to determine the duty cycle of the PWM. I.e.,
  // when the corresponding counter TCNT4 reaches the
  // values in these registers pins OC4A/B will toggle
  // More info: see Atmega32U4 datasheet, section 15.8.2
  TCCR4A |= (1 << PWM4A)  | (1 << PWM4B);
}


void set_fast_pwm_mode(){
  // Enable Fast PWM Mode (WGM41:40 = 0B00)
  TCCR4D |= (0 << WGM41) | (0 << WGM40);
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
