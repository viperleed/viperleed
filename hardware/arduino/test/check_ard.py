# for testing only, delete when done testing
import time
import pyautogui as pya

def main():
    num_relay = 0
    num_LM35 = 0
    count = 1000
    x, y = pya.locateCenterOnScreen('arduino.png', confidence=0.99)
    location = (x, y, 800, 400)
    x, y = pya.locateCenterOnScreen('send.png', confidence=0.99)
    for i in range(count):
        # time.sleep(1)
        pya.click(x, y)
        a = pya.locateOnScreen('relay.png', confidence=0.99, region=location)
        if a is not None:
            num_relay += 1
        b = pya.locateOnScreen('LM35.png', confidence=0.99, region=location)
        if b is not None:
            num_LM35 += 1
    print('relay: ', num_relay)
    print('LM35: ', num_LM35)
    
if __name__ == '__main__':
    main()
