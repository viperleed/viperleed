

#define F_CLK_T4 16e6             // Timer/Counter4 clock = 16 MHz
#define COIL_1 PC7
#define COIL_2 PB6



void setup()
{ 
  pwm_config(20000, 0.625, 0);    
  pwm_config(20000, 0.25, 1);


  pinMode(COIL_1, OUTPUT);             // Define PC7 (OC4A) as output
  pinMode(COIL_2, OUTPUT);             // Define PC7 (OC4A) as output
}

// Register names 'OCRnx' contain device number 'n' (where 'n' is Timer/Counter n) and Output Compare unit 'x' (where 'x' is A/B/C)
// E.g. OCR4A = Output Compare register on Timer/Counter 4, channel A 
// Function call 'pwm_config' makes use of 10-bit Timer/Counter4
// Arguments:
// 'freq': PWM output frequency, range: 15625 Hz to ~ 1 MHz 
// 'duty_cycle': Duty cycle, range 0.0 to 1.0
// 'channel': Timer/Counter4 channel A or B (pins OC4A or OC4B); Values = 0 / 1 (channel A / B)
void pwm_config(double freq, double duty_cycle, uint8_t channel)
{
  // 'uint16_t' is sufficient because TC4H:OCR4C is 10 bits wide
  uint16_t temp;

  // Pins OC4A, OC4B: cleared on compare match (TCNT = OCR4A/B), set when TCNT = 0x000; Enable PWM output channels A and B 
  TCCR4A = (1 << COM4A1) | (0 << COM4A0) | (1 << COM4B1) | (0 << COM4B0) | (1 << PWM4A) | (1 << PWM4B);

  // Timer/Counter4 prescaler = 1; Timer/Counter4 clock frequency = system clock frequency / 1 (f_clk_T4 = 16 MHz)
  TCCR4B = (0 << CS43) | (0 << CS42) | (0 << CS41) | (1 << CS40);

  // Disable Enhanced Compare/PWM mode (ENHC4 = 0)
  TCCR4E = (0 << ENHC4);

  // Enable Fast PWM Mode (WGM41:40 = 0B00) 
  TCCR4D = (0 << WGM41) | (0 << WGM40);   


  temp = F_CLK_T4 / freq;

  // Set PWM frequency to 20 kHz (TC4H:OCR4C = 799 gives 20 kHz exactly; datasheet formula off by one for fast PWM mode)
  TC4H = (temp - 1) >> 8;
  OCR4C = (temp - 1) & 255;

  switch(channel)
  {
    case 0:
      // Set PWM duty cycle for channel A: duty_cycle = TC4H:OCR4A / TC4H:OCR4C  
      TC4H = (uint16_t)((temp - 1) * duty_cycle) >> 8;
      OCR4A = (uint16_t)((temp - 1) * duty_cycle) & 255;                        
      break;

    case 1:
      // Set PWM duty cycle for channel B: duty_cycle = TC4H:OCR4B / TC4H:OCR4C
      TC4H = (uint16_t)((temp - 1) * duty_cycle) >> 8;
      OCR4B = (uint16_t)((temp - 1) * duty_cycle) & 255; 
      break;
  }
}


void loop()
{ 
  while(1)
  {    
    delay(1000);
    
    // Serial.print("Register value of CLKSTA: ");       // Clock Status Register
    // Serial.println(CLKSTA, BIN);

    // Serial.print("Register value of CLKPR: ");        // Clock Prescaler Register
    // Serial.println(CLKPR, BIN);

    // Serial.print("Register value of PLLCSR: ");       // PLL Control and Status Register
    // Serial.println(PLLCSR, BIN);

    // Serial.print("Register value of PLLFRQ: ");       // PLL Frequency Control Register
    // Serial.println(PLLFRQ, BIN);


    //Serial.print("Register value of PORTB: ");
    //Serial.println(PORTB, HEX);
  }
}