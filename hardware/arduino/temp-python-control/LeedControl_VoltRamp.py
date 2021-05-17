"""
Created on Fri Jul 31 13:38:53 2020

@author: Bernhard Mayr
@author: Michele Riva (from 2020-11-09)
@author: Florian Dörr
"""
# Python standard modules
import time
import sys, os
import configparser
import glob
import struct
from zipfile import ZipFile

# NON STANDARD. Will try to get rid of these
import serial        # NON STANDARD - maybe one can use this to replace the serial package: https://github.com/wiseman/arduino-serial
import pandas as pd  # NON STANDARD I will try to get rid of this. It's used only for the export to csv
import serial.tools.list_ports  # NON STANDARD

current_path = os.path.dirname(os.path.abspath(__file__))

if 'Camera_libraries' not in sys.path: 
    sys.path.append(os.path.join(current_path, 'Camera_libraries'))

# ViPErLEED
from camera import Camera

# Configuration-File location:
configfile_location = 'Configuration/LeedControl_config.ini'
config = configparser.ConfigParser()
if not config.read(configfile_location):
    raise AttributeError("Couldn't load configuration file ",
                         configfile_location)

PC_OK = int(config['communication_bytes']['PC_OK'])

# order of most significant byte in int.to_bytes(b, 'little/big')
# int.to_bytes(bytes_lengthg, orderofsignificantByte)
sign = 'big'

# serial port object to Arduino
arduino_port = '' 

#======================================
#TODO: Change names to something more explanatory

def sendToArduino(send_bytes):
    """
    Function sends a byte or bytearray to the arduino and masks
    the message with STARTMARKER + message + ENDMARKER
    
    Parameters
    ----------
    send_bytes: int, float or bytearray
    """
    STARTMARKER = int(config['communication_bytes']['STARTMARKER'])
    ENDMARKER = int(config['communication_bytes']['ENDMARKER'])
    
    send_bytes, bytes_length = encodeHighBytes(send_bytes)
    adjusted_sendbyte = (STARTMARKER.to_bytes(1, sign)
                         + bytes_length.to_bytes(1, sign)
                         + send_bytes + ENDMARKER.to_bytes(1, sign))
    arduino_port.write(adjusted_sendbyte)

def encodeHighBytes(in_bytes):
    """
    Function encodes message before sending it to the arduino. Encodes int 
    to byte or to a bytearray. If resulting byte or bytearray contains a 
    SPECIAL_BYTE, STARTMARKER or an ENDMARKER it exchange this msg_byte to two
    bytes with a leading SPECIAL_BYTE and (msg_byte - SPECIAL_BYTE). Example: 
        message = b'\x00\x01\xff\x02' (0,1,255,2)
        SPECIAL_BYTE = b'\xfd' (253)
        start_byte = b'\xff' (255)
        
        encoded_byte = message[2] - SPECIAL_BYTE = 255 - 253 = 2 = b'\x02'
        
        new_message = message[0:2] + SPECIAL_BYTE + encoded_byte + message[3] 
                    =  b'\x00\x01\xfd\x02\x02 (0,1,253,2,2)
    
    Parameters
    ----------
    in_bytes: byte, bytearray or int
    
    Returns
    ---------
    outstring:  bytearray or byte
    n:          length of bytearray or byte, dtype: int
    """
    SPECIAL_BYTE = int(config['communication_bytes']['SPECIAL_BYTE'])
    
    out_bytes = ""
    out_bytes = out_bytes.encode()
    if isinstance(in_bytes, int): 
        n = 1
        x = in_bytes
        x >>= 8
        while x != 0: 
            x >>= 8
            n += 1
        in_bytes = in_bytes.to_bytes(n, sign)
    else:
        n = len(in_bytes)
    for i in range(n):
        if in_bytes[i] >= SPECIAL_BYTE:
            out_bytes += SPECIAL_BYTE.to_bytes(1, sign)
            out_bytes += (in_bytes[i] - SPECIAL_BYTE).to_bytes(1, sign)
        else:
            out_bytes = out_bytes + in_bytes[i].to_bytes(1, sign)
    return(out_bytes, n)

#======================================    


def receiveFromArduino():
    """
    Function receives messages from Arduino
    
    Return
    ----------
    send_bytes: int, float or bytearray
    """
    STARTMARKER = int(config['communication_bytes']['STARTMARKER'])
    ENDMARKER = int(config['communication_bytes']['ENDMARKER'])
    
    received_value = "z"
    message = ''
    t1 = time.time()
    while arduino_port.inWaiting() == 0: 
        if (time.time() - t1) > 30: 
            arduino_port.close()
            raise IOError("Timeout while Reading Arduino. No message received")
    while ord(received_value) != STARTMARKER: 
        received_value = arduino_port.read()
        message = received_value

    while ord(received_value) != ENDMARKER:
        received_value = arduino_port.read()
        message += received_value
    
    return decodeHighBytes(message)
    
    
def decodeHighBytes(in_bytes):
    """
    Function decodes message. If there are SPECIAL_BYTEs in the message, it 
    decodes them and the byte after back to a single byte. Example: 
        message = b'\x00\x01\xfd\x02\x02' (0,1,253,2,2)
        SPECIAL_BYTE = b'\xfd' (253)
        start_byte = b'\xff' (255)
        
        decoded_byte = message[2] + message[3] = 253 + 2 = 255 = b'\xff'
        
        new_message = message[0:2] + decoded_bytes + message[3] 
                    =  b'\x00\x01\xff\x02'
    
    Parameter
    ----------
    in_bytes: bytearray or byte
    
    Return
    ----------
    byte or bytearray
    """
    SPECIAL_BYTE = int(config['communication_bytes']['SPECIAL_BYTE'])
    
    out_bytes = in_bytes[0].to_bytes(1, sign)
    i = 1
    
    if in_bytes[1] == 0: 
        arduino_port.close()
        raise IOError("Arduino sended an ERROR: %s" % 
                      str(in_bytes[2:-3])[2:-1])
    while i < len(in_bytes):
        if in_bytes[i] == SPECIAL_BYTE:
            out_bytes += (in_bytes[i] + in_bytes[i+1]).to_bytes(1, sign)
            i += 1
        else: 
            out_bytes += in_bytes[i].to_bytes(1, sign)
        i += 1
     
    if out_bytes[1] < 4: 
        return int.from_bytes(out_bytes[2:(2+out_bytes[1])], sign)
    elif out_bytes[1] == 4: 
        return struct.unpack('f', out_bytes[2:6])[0]
    
#======================================  

def connectToArduino(port):
    """
    Function connects to Arduino, through given port ID
    The Arduino will return the hardware it can detect to
    the PC
    
    Parameter
    ----------
    port: string
    
    Return
    ----------
    arduino_port: class object serial
    """
    global arduino_port
    PC_CONFIGURATION = int(config['communication_bytes']['PC_CONFIGURATION'])
    FIRMWARE_VERSION_MAJOR = int(config['firmware_version']['FIRMWARE_VERSION_MAJOR'])
    FIRMWARE_VERSION_MINOR = int(config['firmware_version']['FIRMWARE_VERSION_MINOR'])
    hardware_bits = int(config['hardware_bits'])
    
    arduino_port = serial.Serial(port, 115200, timeout = 1)
    if arduino_port.inWaiting() == 0:
        sendToArduino(PC_CONFIGURATION)
        time.sleep(1) #time to settle connection
        hardwareConfiguration = receiveFromArduino()
        if hardwareConfiguration[1] != FIRMWARE_VERSION_MAJOR or hardwareConfiguration[2] != FIRMWARE_VERSION_MINOR:
            print("Versions do not match.")
        if hardwareConfiguration[3] or hardwareConfiguration[4]:
            print("Hardware detected")
            for key, value in hardware_bits.items():
                if value & hardwareConfiguration:
                    print (key)
            return arduino_port
        elif hardwareConfiguration is None:
            print("No hardware configuration returned")
        elif hardwareConfiguration[3] == 0 and hardwareConfiguration[4] == 0:
            print("No hardware detected")
            return arduino_port
    arduino_port.close()
    raise IOError("\nSomething went wrong! \nConnection was successful but no" 
                  " respond from the device! \nT H E  E N D")

#======================================


def initialiseADC(arduino_config):
    """
    Function initialise adc and sends adc configuration settings, defined in
    'ReadArdo_config.ini'
    
    Parameter
    ----------
    arduino_config: dict
    """
    PC_SET_UP_ADCS = int(config['communication_bytes']['PC_SET_UP_ADCS'])
    sendToArduino(PC_SET_UP_ADCS)
    message = []
    measurement_n = int(arduino_config['measurement_counter']).to_bytes(
        2, sign)
    message.append(measurement_n[0])
    message.append(measurement_n[1])
    for key in ('adc0_channel', 'adc1_channel'): 
        message.append(int(arduino_config[key]))
    sendToArduino(message) 
    #print("Gain in ADC initialisation is:", receiveFromArduino())
    if receiveFromArduino() == PC_OK: 
        print("ADC Initalisation: DONE!")
        print("-----------------------")
    else:
        identifyError()
        arduino_port.close()
        raise IOError("ADC Initalisation: FAILED!")

#======================================


def serial_ports():
    """ 
    Picks Arduino connection from all serial connections

    Raises
    -------
    EnvironmentError: On unsupported or unknown platforms
    
    Returns
    -------
    result: string
    """
    results = []
    
    if sys.platform.startswith('win'):
        ports = list(serial.tools.list_ports.comports())
        for port in ports: 
            if 'Arduino Micro' in port[1]: results.append(port[0])
    elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # this excludes your current terminal "/dev/tty"
        ports = glob.glob('/dev/tty[A-Za-z]*')
    elif sys.platform.startswith('darwin'):
        ports = glob.glob('/dev/tty[A-Za-z]*')
    else:
        raise EnvironmentError('Unsupported platform')

    result_port = ''
    
    for result in results:
        try:
            s = serial.Serial(result)
            s.close()
            result_port = result
        except (OSError, serial.SerialException):
            pass
        except ValueError: 
            raise IOError("Couldn't find port to Arduino!")
    if not results: 
        raise IOError("Couldn't find port to Arduino!")
    if not result_port:
        raise IOError("Unable to connect to Arduino, try disconnecting and "
                      "connecting again.")
    return result_port

#======================================
def setVoltageAndMeasure(energy, settle_time = 0):
    """
    Function initialise dac and sends dac value to arduino. 
    
    Parameter
    ----------
    energy: float, int
    settle_time: int
    """
    v_ref_dac = float(config['measurement_settings']['v_ref_dac'] )
    PC_SET_VOLTAGE = int(config['communication_bytes']['PC_SET_VOLTAGE'])
    dac_value_temp = int(round(energy * 65536/(v_ref_dac*2)))
    if(dac_value_temp >= 65536): 
        dac_value_temp = 65535
    dac_send = dac_value_temp.to_bytes(2,sign) + settle_time.to_bytes(2,sign)
    # print("DAC_send = " ,dac_value_temp)
    # print("DAC_tobytes =" ,dac_value_temp.to_bytes(2, sign), " t= ",round(1000*time.time()))
    sendToArduino(PC_SET_VOLTAGE)
    sendToArduino(dac_send)
    # dac_send = sendToArduino(dac_send)
    if receiveFromArduino() == PC_OK:
        print("Voltage is set, triggering")
    else:
        identifyError()
        arduino_port.close()
#======================================

def initialiseAutoGain():
    """
    Function starts autogain-function on arduino. Autogain-funtion measures 
    adc value and sets gain, that measured value is bigger than 25% of the 
    measurment range. Starts with gain 1 (no gain). Maximum gain is 128
    """
    PC_AUTOGAIN = int(config['communication_bytes']['PC_AUTOGAIN'])
    sendToArduino(PC_AUTOGAIN)
    print("Autogain: IN PROCESS...")
    
    if receiveFromArduino() == PC_OK:
        # print("Gain found =",pow(2,int(receiveFromArduino())))#for debugging, can be deleted
        print("Gain found =", int(receiveFromArduino()))
    else:
        identifyError()
        arduino_port.close()
        raise IOError("Arduino doesnt react when python tries to initialise "
                      "the Autogain Function, Program stops")

#======================================

def createCSV(multilist, path):
    """
    Creates csv with columns:
        Index | ADC_value [V] | DAC_value [V] | Time [sec]
    safes csv in path
    """
    if not os.path.exists(path):
        os.makedirs(path)
        print("New directory with new file was created")
    i = 0
    while os.path.exists(os.path.join(path, "sample%s.csv" % i)):
        i += 1
    #timestamp = datetime.now()
    df = pd.DataFrame(multilist, columns=['ADC_VAL', 'DAC_Val', 'DeltaTime'])
    print("-----------------------") 
    print(df)
    print("-----------------------") 
    df.to_csv(os.path.join(path, "sample%s.csv" % i), index=True)
    print("CSV File was created")

#======================================

def TUI():
    """
    Little terminal user interface for defining the Startenergy, Endenergy 
    and the DeltaEnergy. Main loop generates a ramp from this values and sends 
    values to DAC. 
    
    Return
    ------
    startenergy, endenergy, deltaenergy: float
    """
    while True:
        max_energy = 2*float(config['measurement_settings']['v_ref_dac'])
        start_energy = -1
        while start_energy < 0 or start_energy > max_energy:
            start_energy = float(input("\nPlease enter the energy to START "
                                      "with in Volt. Allowed are values between "
                                      f"0 V and {(max_energy)*100:.1f} V:\n"))/100
            if start_energy < 0 or start_energy > max_energy:
                print("Invalid")
        end_energy = start_energy - 1
        while end_energy < start_energy or end_energy > max_energy:
            end_energy = float(
                input("Please enter the Energy to END with in Volt:\n"))/100
            if end_energy < start_energy or end_energy > max_energy:
                print("Invalid")
        delta_energy = -1
        while delta_energy < 0:
            delta_energy =  float(
                input("Please enter the DELTA energy in Volt:\n"))/100
            if delta_energy < 0:
                print("Invalid")
        print("==================================")
        print(f"Start Energy: {100*start_energy:.1f} V")
        print(f"End Energy: {100*end_energy:.1f} V")
        print(f"Delta Energy: {100*delta_energy:.2f} V")
        print("==================================")
        choose = input("Are you satisfied with your desired values [Y/n]?")
        if choose.lower() == 'y':
            print("Values approved, start with test NOW!")
            print("==================================")
            return start_energy, end_energy, delta_energy 
        print("Values not approved, again...")


def pack_measurements(zip_path):  # MR: NEEDS TO BE COMPLETED
    """
    """
    # Prepare files to pack together, already taking into account their absolute
    # path

    # 1) configuration file
    fnames = [os.path.join(current_path, configfile_location)]

    # 2) csv file
    fnames.extend(f for f in os.listdir(config['measurement_settings']['path'])
                  if ".csv" in f)

    # 3) image files
    fnames.extend(f for f in os.listdir(config['camera_settings']['image_path'])
                  if ".tif" in f)
    
    # TEMP: will use a simple incremental number for a moment for the naming of
    # the zip files
    i = 0
    while os.path.exists(os.path.join(path, str(i), ".zip")):
        i += 1
    
    with ZipFile(os.path.join(zip_path, str(i), ".zip")) as zip_file:
        for f in fnames:
            zip_file.write(f, os.path.basename(f))
            
def adcCalibration(arduino_config):
    """
    Send PC_CALIBRATION to arduino which does the calibration for all gains for the
    selected channels and saves the values for later use.
    
    Parameter
    ----------
    arduino_config: dict
    """
    PC_CALIBRATION = int(config['communication_bytes']['PC_CALIBRATION'])
    sendToArduino(PC_CALIBRATION)
    message = []                                               
    for key in ('update_rate', 'adc0_channel', 'adc1_channel'): 
        message.append(int(arduino_config[key]))
    sendToArduino(message) 
    if receiveFromArduino() == PC_OK: 
        print("ADC calibration: DONE!")
        print("-----------------------")
    else:
        identifyError()
        arduino_port.close()
        raise IOError("ADC calibration: FAILED!")


def identifyError():
    """
    Is called upon when PC_ERROR is returned from the Arduino to the PC. It compares 
    the trace back byte to the ones saved in the config and prints the according key
    """
    arduino_states = int(config['arduino_states'])
    error_bytes = int(config['error_bytes'])
    ErrorMessage = receiveFromArduino()
    for key, value in arduino_states.items():
        if value == ErrorMessage[1]:
            print (key)
    for key, value in error_bytes.items():
        if value == ErrorMessage[2]:
            print (key)

def main(energies=None):
    global arduino_port
    
    PC_RESET = int(config['communication_bytes']['PC_RESET'])
    
    adc_autogain_config = config['adc_autogain_configuration']
    adc_start_config = config['adc_start_configuration']
    v_ref_adc = float(config['measurement_settings']['v_ref_adc'])
    dac_first_settletime = int(config[
        'measurement_settings']['dac_first_settle_time'])
    dac_settletime =  int(
        config['measurement_settings']['dac_settle_time'])
    path = str(config['measurement_settings']['path'])
    dac_value_end = float(config['measurement_settings']['dac_value_end'])
    camera_settle_time = int(config['camera_settings']['camera_settle_time'])
    live_mode = (config['camera_settings']['live_mode'] == 'True')
    I0_conversion_factor = int(config['leed_hardware']['I0_conversion_factor'])
    
    start_energy, end_energy, delta_energy = TUI()
    start_energy = 0.5
    delta_energy = 0.005#1e-16#
    end_energy = 4.9#0.5 + 881*delta_energy*3#
    
    #Connect with Camera and initialize it
    CameraObject = Camera.camera_from_ini(
        config['camera_settings']['class_name'])
    CameraObject.initialize(config['camera_settings'])

    # Connect with Arduino, reset it and flush input buffer
    portname = serial_ports()
    arduino_port = connectToArduino(portname)
    arduino_port.flushInput()

    # Calibrate and initialize ADCs, DAC and auto-gain,
    # then re-initialize the ADC with the right gain
    adcCalibration(adc_start_config)
    initialiseADC(adc_autogain_config)
    setVoltageAndMeasure(start_energy, dac_first_settletime)
    initialiseAutoGain()
    actual_energy = start_energy
    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    dac_value_csv = []
    timestamp = []
    timestart = time.time()
    
    # This keyboard interrupt thing will not exist. Process interruption will
    # be handled with a signal
    # TODO: I have no clue what this does but I removed the call for a measurement as setVoltageAndMeasure already does that. We need to redo this.
    try: 
        while actual_energy <= end_energy:
            print("Energy:%.2f V" % actual_energy)
            # if camera is used in "exposure" mode
            while not live_mode: 
                if not CameraObject.callback_data.dac_busy:
                    setVoltageAndMeasure(actual_energy, dac_settletime)
                    break
            if not live_mode:
                time.sleep(camera_settle_time*1e-3)
                CameraObject.callback_data.process_params[
                    'filename_energy'] = "%.3f" % round(actual_energy, 3)
                # print(actual_energy)
                t1 = time.time()
                CameraObject.trigger_now()
                print("Time:", (time.time()-t1))
            else: 
                pass
                # CameraObject.snap_image()
            #print('Gain: %i \nDAC-Value: %f' % (int(receiveFromArduino()),
            #                                   actual_energy))
            print('nDAC-Value: %f' % (actual_energy))
            adc0_value = receiveFromArduino()
            adc1_value = receiveFromArduino()
            adc2_value = receiveFromArduino()
            #adc0 is ADC0, adc1 is ADC1 and adc2 is the LM35 in the Arduino code
            print("ADC0_VAlUE:", adc0_value*I0_conversion_factor)
            print("ADC1_VAlUE:", adc1_value)
            print("ADC2_VAlUE:", adc2_value,
                  " t = ", round(1000*time.time()))
            print("#########################")
            adc0_value_csv.append(adc0_value*I0_conversion_factor)
            adc1_value_csv.append(adc1_value)
            adc2_value_csv.append(adc2_value)
            dac_value_csv.append(actual_energy)
            time_temp = time.time() - timestart
            timestamp.append(time_temp)
            actual_energy += delta_energy

    # finally: 
    except KeyboardInterrupt:
        CameraObject.stop_camera()
        #sendToArduino(PC_RESET)
        #Reset would clear calibration: not necessary
        setVoltageAndMeasure(dac_value_end) 
        arduino_port.close()
        
    time.sleep(3)  # To process the last frame, the callback needs some time
    CameraObject.stop_camera()
    setVoltageAndMeasure(dac_value_end)
    createCSV(list(zip(adc0_value_csv, adc1_value_csv, adc2_value_csv, dac_value_csv, timestamp)), path) 
    arduino_port.close()
    print('Arduino Disconnected')
    
    # before terminating, pack all the necessary data into a single zip file
    # pack_measurements()

def do_nothing_on_purpose():  # MR: just a mock I needed for my presentation
    import numpy as np
    
    start_energy, end_energy, delta_energy = TUI()
    actual_energy = start_energy
    
    time.sleep(1) #time to settle connection
    print("Arduino Initialisation: DONE!")
    
    time.sleep(0.5) #time to settle connection
    print("ADC Initalisation: DONE!")
    print("-----------------------")
    time.sleep(0.2) #time to settle connection
    print("Autogain: IN PROCESS...")
    
    time.sleep(0.2) #time to settle connection
    print("ADC Initalisation: DONE!")
    print("-----------------------")
    
    while actual_energy <= end_energy/3:
        print("Energy:%.2f V" % (actual_energy*100))
        t1 = time.time()
        time.sleep(0.1)
        print("Time:", (time.time()-t1))
        print('Gain: %i \nDAC-Value: %f' % (64, actual_energy*100))
        print("ADC_VAlUE:", np.random.rand(1) * 0.49,
                      " t = ", round(1000*time.time()))
        print("#########################")
        actual_energy += delta_energy
    
    time.sleep(3)
    print('Arduino Disconnected')



if __name__ == '__main__':
    t1 = time.time()
    for i in range(1):
        main()
        print('#########################')
        print('DONE!')
        print("Time for whole procedure:", (time.time() - t1))
