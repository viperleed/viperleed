/*
ViPErLEED - Driver for DRV8842.
---------------------
Authors: Stefan Mitterhöfer
26.11.2024, SM: First version
---------------------
*/

#ifndef _VIPERLEED_DRV8842
#define _VIPERLEED_DRV8842

class DRV8842 : public MotorDriver {
    public:
        DRV8842(byte enable_pin)
        : MotorDriver(enable_pin) {}

        DRV8842(byte chip_select_pin, byte enable_pin)
        : MotorDriver(chip_select_pin, enable_pin) {}

//        void setup() override {
//            pinMode(_ENABLE_PIN, OUTPUT);     // Use for DRV8842 enable pin 'nRESET'
//            digitalWrite(_ENABLE_PIN, HIGH);  // 'nRESET' == 1: Activate DRV8842 output drivers
//        }

    private:


};

#endif  // _VIPERLEED_DRV8842
