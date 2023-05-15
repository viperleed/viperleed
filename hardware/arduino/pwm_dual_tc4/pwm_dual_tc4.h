/*
Test firmware for Arduino dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/


#ifndef _VIPERLEED_B_FIELD_COMP
#define _VIPERLEED_B_FIELD_COMP


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
                                
                                
                                
#endif
