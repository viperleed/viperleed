#include <SPI.h>

#define TLE_CHIPSELECT 11 // CS line on IO11/PB7


// PWM period = 50 µs / PWM frequency = 20 kHz 
void setup()
{ 
  // Enable Fast PWM Mode (PWM4X = 1, WGM41..40 = 0), TOP = OCR4C
  TCCR4D &= ~ 0x3;

  TCCR4A = (1 << COM4A0) + (1 << COM4B0) + (1 << PWM4A) + (1 << PWM4B);           // Toggle OC4A (pin PC7) and OC4B (PB6) at compare match
  TCCR4B = (0 << CS43) + (0 << CS42) + (1 << CS41) + (1 << CS40);                 // prescaler = 4; 

  OCR4C = 200;                    // PWM period = T_clk * prescaler * OCR4C;  T_clk = 62.5 ns (Arduino Micro)       

  OCR4A = 100;                    // OC4A PWM duty cycle = OCR4A / OCR4C
  OCR4B = 100;                    // OC4B PWM duty cycle = OCR4B / OCR4C

  DDRC |= (1 << PC7);             // Define PC7 as output (Pin PC7 = IO13)
  DDRB |= (1 << PB6);             // Define PB6 as output (Pin PB6 = IO10)

  Serial.begin(9600);
  //Serial.setTimeout(100);

  pinMode(TLE_CHIPSELECT, OUTPUT);  

  SPI.begin();                    // SPI.begin initializes the SPI bus (defines SCK and MOSI as output pins)
                                  // Alternatively, use pinMode() on SCK and MOSI below
  //pinMode(SCK, OUTPUT);           // SCK = pin PB1
  //pinMode(MOSI, OUTPUT);          // MOSI = pin PB2
  //pinMode(MISO, INPUT);           // MISO = pin PB3
} 


// uint8_t data[] = {0x0000};

void loop()
{ 
  // while(1)
  // {
  //   Serial.println("\nType 'go' to read register DDRB: ");
  //   while (Serial.available() == 0) {}     //wait for data available

  //   String teststr = Serial.readString();  //read until timeout
  //   teststr.trim();                        // remove any \r \n whitespace at the end of the String
  //   if (teststr == "go") 
  //   {
  //     Serial.print("Register value of DDRB: ");
  //     Serial.println(DDRB, HEX);
  //   } else 
  //   {
  //     Serial.println("Something else");
  //   }
  // }


  while(1)
  {
    delayMicroseconds(200);
    //Serial.print("Serial debug\n\n");    

    // Note: The TLE7209 always operates in slave mode (MOSI = SDI, MISO = SDO)
    // SPI Mode 1: CPOL = 0, CPHA = 1
    digitalWrite(TLE_CHIPSELECT, LOW);
    delayMicroseconds(5);
    SPI.beginTransaction(SPISettings(1000000, MSBFIRST, SPI_MODE1));

    SPI.transfer(0x00);
    SPI.transfer(0x00);

    digitalWrite(TLE_CHIPSELECT, HIGH);
    SPI.endTransaction();  

    delayMicroseconds(10);

  }
}