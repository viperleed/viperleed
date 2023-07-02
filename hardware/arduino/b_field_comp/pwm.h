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
   

// 'pwm_period': How many Timer/Counter4 clock ticks until TCNT4 rolls over.
// This essentially sets the frequency of the PWM. 'uint16_t' is sufficient 
// because TC4H:OCR4C is 10 bits wide. Initialized inside 'set_pwm_frequency'.
uint16_t pwm_period;

#define POSITIVE_CURRENT  1
#define NEGATIVE_CURRENT -1

                                // register OCR4C and from the choice of a
                                // value of 1 in set_pwm_clock_prescaler
#define F_CPU_CLK 16e6                                    // Arduino Micro CPU clock = 16 MHz
#define F_PWM_MIN (F_CPU_CLK / (16384 * 1024.0))          // F_PWM_MIN = 0.954 Hz
#define F_PWM_MAX (F_CPU_CLK / (1 * 4.0))                 // F_PWM_MAX = 4 MHz
enum TC4_PWM_CHANNEL { TC4_PWM_CH_A, TC4_PWM_CH_B, TC4_PWM_CH_D };

const uint16_t TC4_CLK_PRESCALER[] = { 0, 1, 2, 4, 8, 16, 32, 64, 128, 256,
                                       512, 1024, 2048, 4096, 8192, 16384};

#endif
