/*
Test firmware for Arduino dual channel PWM generation
---------------------
Author: Michele Riva, Christoph Pfungen
Date: 15.05.2023
---------------------
*/


#ifndef _VIPERLEED_B_FIELD_PWM
#define _VIPERLEED_B_FIELD_PWM

#define FAST_PWM_CH_1_REG 0xD2      // Equals the address of register OCR4D, see sec. 31 Register Summary
#define FAST_PWM_CH_2_REG 0xD0      // Equals the address of register OCR4B, see sec. 31 Register Summary
   

// pwm_clock_divider: how many intervals the counters count before
// rolling over. This essentially sets the frequency of the PWM.
// 'uint16_t' is sufficient because TC4H:OCR4C is 10 bits wide.
uint16_t pwm_clock_divider;  // Use set_pwm_frequency for setting!

#define POSITIVE_CURRENT  1
#define NEGATIVE_CURRENT -1

#define F_CLK_T4     16000      // Timer/Counter4 clock = 16 MHz
#define PWM_MIN_FREQ 15.625     // Value comes from maximum value for
                                // register OCR4C and from the choice of a
                                // value of 1 in set_pwm_clock_prescaler

#endif
