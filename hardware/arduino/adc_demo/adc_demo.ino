

#include <Arduino.h>
#include <LiquidCrystal.h>
#include "adc_demo.h"

#define DEBUG true
#define ADC_INPUT_0 A0        // Pin PF7 on the Arduino Micro board


const int rs = 12, en = 11, d4 = 5, d5 = 4, d6 = 3, d7 = 2;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);

void setup() 
{
  // Use AVcc as reference voltage, result is right-adjusted, single-ended input on pin ADC0
    ADMUX = (0 << REFS1) | (1 << REFS0) | (0 << ADLAR) | (0 << MUX4) 
          | (0 << MUX3)  | (0 << MUX2)  | (0 << MUX1)  | (0 << MUX0);

  // ADC prescaler value = 128
  ADCSRA = (0 << ADEN) | (0 << ADSC) | (0 << ADATE) | (0 << ADIF) 
         | (0 << ADIE)  | (1 << ADPS2)  | (1 << ADPS1)  | (1 << ADPS0);    // TODO: apply correct ADC clock prescaler

  // Set up ADC in free running mode / manual trigger
  ADCSRB = (0 << ADHSM) | (0 << ACME) | (0 << MUX5) 
         | (0 << ADTS3)  | (0 << ADTS2)  | (0 << ADTS1)  | (0 << ADTS0);

  // Disable digital input on pin ADC0
  DIDR0 = 1 << ADC0D;
  DIDR2 = 0;

  // Enable ADC
  ADCSRA |= 1 << ADEN;


  pinMode(ADC_INPUT_0, INPUT);
  //pinMode(PWM_OUTPUT, OUTPUT);

  // Set up the LCD's number of columns and rows:
  lcd.begin(16, 2);

  #if DEBUG
    Serial.setTimeout(100);
    Serial.begin(9600);           // Open serial port, set data rate to 9600 bps
  #endif
}

void loop() 
{
  delay(500);

  ADCSRA |= 1 << ADSC;            // Start A/D conversion
  while(bitRead(ADCSRA, ADSC));   // Wait while conversion is in progress

  lcd.clear();
  lcd.setCursor(0, 0);
  lcd.print(ADC);
}
