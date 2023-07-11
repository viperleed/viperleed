/**Useful functions for Arduino, not included in STDlib

@author: Michele Riva
@author: Michael Schmid
@author: Florian Doerr
@author: Christoph Pfungen
**/

#ifndef _VIPERLEED_ARDUINO_UTILS
#define _VIPERLEED_ARDUINO_UTILS


#include <Arduino.h>  // For pin definitions and macros


void setChipSelectHigh(byte ioPin) {
    /**
    Set a digital output used as a chip select signal from
    the default high-impedance state to high (=unselected),
    without a glitch to the low state.
    **/
    pinMode(ioPin, INPUT_PULLUP);
    digitalWrite(ioPin, HIGH);
    pinMode(ioPin, OUTPUT);
}


// If pin definitions on the ATmega32U4 should ever change, 
// the assignments below would have to change accordingly.
uint8_t pin_to_tc4_reg_addr(uint8_t pwm_pin) {
    switch(pwm_pin) {
      case  6: return 0xD2;
      case 10: return 0xD0;
      case 13: return 0xCF;
    }
}


// If pin definitions on the ATmega32U4 should ever change, 
// the assignments below would have to change accordingly.
uint8_t pin_to_tc4_channel(uint8_t pwm_pin) {
    switch(pwm_pin) {
      case  6: return TC4_PWM_CH_D;
      case 10: return TC4_PWM_CH_B;
      case 13: return TC4_PWM_CH_A;
    }
}


uint16_t bigger16(uint16_t a, uint16_t b) {
  return (a > b) ? a : b;
}

uint16_t biggest16(uint16_t a, uint16_t b, uint16_t c) {
  return bigger16(a, bigger16(b, c));
}

/** Gets the median of three numbers */
uint16_t getMedian16(uint16_t a0, uint16_t a1, uint16_t a2) {
  uint16_t maximum = biggest16(a0, a1, a2);
  if (maximum == a0) return bigger16(a1, a2);
  if (maximum == a1) return bigger16(a0, a2);
  else return bigger16(a0, a1);
}

uint16_t analogReadMedian(byte pin) {
  /**
  Read the Arduino ADC for a given pin
  (A0...) three times and return the median
  */
  uint16_t a0 = analogRead(pin);
  uint16_t a1 = analogRead(pin);
  uint16_t a2 = analogRead(pin);
  return getMedian16(a0, a1, a2);
}

#endif
