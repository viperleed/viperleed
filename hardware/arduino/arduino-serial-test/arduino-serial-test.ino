/*Simple serial-line tester sketch for Arduino Micro.

Reads from serial line and echoes back the message after
a full one was received.

Created: 2021-06-30
Author: Michele Riva
*/

#include <Arduino.h>

#define MSG_START  0x00
#define MSG_END    0xff

#define MONITOR false

byte serialBuffer[16];  // Serial buffer
byte numBytesRead;
bool receiving;


void setup(){
    // Do nothing, as the Serial module does not need a begin on the Micro
    #if MONITOR
        Serial.begin(115200);
        Serial.println("--- Start Serial Monitor SEND_RCVE ---");
        Serial.println("(Decimal)(Hex)(Character)");
        Serial.println();
    #endif
}

void loop(){
    if (not receiving){
        delayMicroseconds(100);
    }

    readAndEchoSerial();
}

void readAndEchoSerial(){
    /*Read from Serial and echo back after a full message is read*/
    if (not Serial.available()) return;

    while (Serial.available()){
        byte character = Serial.read();

        #if MONITOR
            Serial.print(character);
            Serial.print("        ");
            Serial.print(character, HEX);
            Serial.print("       ");
            Serial.print(char(character));
            Serial.println();
        #endif

        if (character == MSG_START){
            numBytesRead = 0;
            receiving = true;
        }

        if (receiving){
            serialBuffer[numBytesRead] = character;
            numBytesRead++;
        }

        if (character == MSG_END){
            #if MONITOR
                Serial.print("E");
                Serial.println();
            #endif
            receiving = false;
            echoSerial();
            return;
        }

        delayMicroseconds(100);
    }
}

void echoSerial(){
    #if MONITOR
      Serial.print("Echoing: ");
    #endif
    for (int i; i < numBytesRead; i++)
        Serial.write(serialBuffer[i]);
    #if MONITOR
      Serial.println();
    #endif
}
