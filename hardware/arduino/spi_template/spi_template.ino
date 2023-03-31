#include <SPI.h>

#define TLE_CHIPSELECT 11 // CS line on IO11/PB7


// PWM period = 50 µs / PWM frequency = 20 kHz 
void setup()
{ 
 
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