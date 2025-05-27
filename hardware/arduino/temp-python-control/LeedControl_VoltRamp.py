"""
Created on Fri Jul 31 13:38:53 2020

@author: Bernhard Mayr
@author: Michele Riva (from 2020-11-09)
@author: Florian Dörr (from 2021-04-26)
"""
# Python standard modules
import time
import sys
import os
import configparser
import glob
import struct
import ast
from zipfile import ZipFile
import matplotlib.pyplot as plt

import numpy as np
from numpy.polynomial.polynomial import Polynomial
# import scipy
from scipy.fft import rfft, rfftfreq
from PyQt5 import (QtSerialPort as qts,
                   QtCore as qtc,
                   QtWidgets as qtw)
# NON STANDARD. Will try to get rid of these
import pandas as pd  # NON STANDARD I will try to get rid of this. It's used only for the export to csv

current_path = os.path.dirname(os.path.abspath(__file__))
dll_path = os.path.join(current_path, 'Camera_libraries')
camera_path = dll_path

try:
    os.add_dll_directory(dll_path)
except AttributeError:
    pass
if camera_path not in sys.path:
    sys.path.append(camera_path)

# ViPErLEED
from camera import Camera

# Configuration-File location:
configfile_location = 'Configuration/LeedControl_config.ini'
config = configparser.ConfigParser()
if not config.read(configfile_location):
    raise AttributeError("Couldn't load configuration file ",
                         configfile_location)

# order of most significant byte in int.to_bytes(b, 'little/big')
# int.to_bytes(bytes_lengthg, orderofsignificantByte)
sign = 'big'

PC_OK = config.getint('communication_bytes', 'PC_OK').to_bytes(1, sign)

# QSerialPort object to Arduino
arduino_port = qts.QSerialPort()

# ======================================

def send_to_arduino(send_bytes):
    """
    Function sends a byte or bytearray to the arduino and masks
    the message with STARTMARKER + message + ENDMARKER

    Parameters
    ----------
    send_bytes: int, float, byte or bytearray
    """
    STARTMARKER = config.getint('communication_bytes', 'STARTMARKER')
    ENDMARKER = config.getint('communication_bytes', 'ENDMARKER')

    send_bytes, bytes_length = encode_high_bytes(send_bytes)
    adjusted_sendbyte = (STARTMARKER.to_bytes(1, sign)
                         + bytes_length.to_bytes(1, sign)
                         + send_bytes + ENDMARKER.to_bytes(1, sign))
    print("Sending: ", adjusted_sendbyte, "with msg_length =", bytes_length)
    if arduino_port.write(adjusted_sendbyte) < 0:
        raise RuntimeError("Serial.write did not succeed")

def encode_high_bytes(in_bytes):
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
    SPECIAL_BYTE = config.getint('communication_bytes', 'SPECIAL_BYTE')

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
    return out_bytes, n

def bytes_to_float(bytes_in):
    """Convert four bytes to float.

    Parameters
    ----------
    bytes_in : bytes or bytearray
        Should have length 4. It is assumed that the bytes are in
        big-endian format.

    Returns
    -------
    float
    """
    if (not hasattr(bytes_in, '__len__')
            or not isinstance(bytes_in, (bytes, bytearray))):
        raise TypeError("bytes_to_float: invalid type "
                        f"{type(bytes_in).__name__}. Expected 'bytes',"
                        "or 'bytearray'")
    if len(bytes_in) != 4:
        raise ValueError("bytes_to_float: Invalid number of bytes "
                         f"({len(bytes_in)}). Expected 4.")
    return struct.unpack(">f", bytes_in)[0]

# ======================================
def adc_measure_only():
    """
    Send  PC_MEASURE_ONLY to arduino and initialise measurement of adc
    Not used in our current code
    """
    PC_MEASURE_ONLY = config.getint('communication_bytes', 'PC_MEASURE_ONLY')
    send_to_arduino(PC_MEASURE_ONLY)
    if receive_from_arduino() == PC_OK:
        print("Measurements starting")
    else:
        identify_error()
        arduino_port.close()
        raise IOError("Arduino does not reach the number of measurements to do")
# ======================================

def receive_from_arduino():
    """
    Function receives messages from Arduino

    Return
    ----------
    send_bytes: int, float or bytearray
    """
    STARTMARKER = config.getint('communication_bytes', 'STARTMARKER')
    ENDMARKER = config.getint('communication_bytes', 'ENDMARKER')
    arduino_timeout = config.getfloat('measurement_settings', 'arduino_timeout')
    received_value = "z"
    message = ''
    t1 = time.time()
    while not arduino_port.bytesAvailable():
        if (time.time() - t1) > arduino_timeout:
            arduino_port.close()
            raise IOError("Timeout while Reading Arduino. No message received")
    while ord(received_value) != STARTMARKER:
        received_value = arduino_port.readData()
        message = received_value

    while ord(received_value) != ENDMARKER:
        received_value = arduino_port.read()
        message += received_value

    # TODO: check that the message length is appropriate

    # print("Got message: ", repr(message))
    decoded = decode_high_bytes(message)
    # print("Decoded: ", repr(decoded))

    return decoded


def decode_high_bytes(in_bytes):
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
    SPECIAL_BYTE = config.getint('communication_bytes', 'SPECIAL_BYTE')

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
    return out_bytes[2:-1]

# ======================================

def connect_to_arduino(port):
    """
    Function connects to Arduino, through given port ID
    The Arduino will return the hardware it can detect to
    the PC

    Parameter
    ----------
    port: string

    Return
    ----------
    arduino_port: class object QSerialPort
    """
    global arduino_port
    PC_CONFIGURATION = config.getint('communication_bytes', 'PC_CONFIGURATION')
    FIRMWARE_VERSION = config['firmware_version']['FIRMWARE_VERSION']
    hardware_bits = config['hardware_bits']

    # arduino_port = serial.serial(port, 115200, timeout = 1)
    arduino_port.setPort(port)
    arduino_port.open(arduino_port.ReadWrite)
    arduino_port.setBaudRate(int(400e6))
    arduino_port.setDataTerminalReady(True)
    # if arduino_port.inWaiting() == 0:  # old solution for PySerial
    # if arduino_port.waitForReadyRead():
    if not arduino_port.bytesAvailable():
        send_to_arduino(PC_CONFIGURATION)
        time.sleep(1)  # time to settle connection
        hardware_config = receive_from_arduino()
        hardwareVersion = f"{hardware_config[0]}.{hardware_config[1]}"
        if hardwareVersion != FIRMWARE_VERSION:
            print(f"Versions do not match: got {hardwareVersion}",
                  "expected", FIRMWARE_VERSION)
        else:
            print(f"Correct version found: v{hardwareVersion}")
        if any(hardware_config[2:]):
            print("Hardware detected")
            for key, value in hardware_bits.items():
                if int(value) & int.from_bytes(hardware_config[2:], sign):
                    print (key)
            return arduino_port
        elif hardware_config is None:
            print("No hardware configuration returned")
        elif hardware_config[2] == 0 and hardware_config[3] == 0:
            print("No hardware detected")
            return arduino_port
    arduino_port.close()
    raise IOError("\nSomething went wrong! \nConnection was successful but no"
                  " response from the device! \n")

# ======================================


def initialise_adcs(arduino_config):
    """
    Function initialise adcs and sends adc configuration settings, defined in
    'LeedControl_config.ini'

    Parameter
    ----------
    arduino_config: dict
    """
    PC_SET_UP_ADCS = config.getint('communication_bytes', 'PC_SET_UP_ADCS')
    send_to_arduino(PC_SET_UP_ADCS)
    message = []
    measurement_n = int(arduino_config['measurement_counter']).to_bytes(2, sign)
    # print(measurement_n) prints b'\x00\n'
    message.append(measurement_n[0])
    message.append(measurement_n[1])
    for key in ('adc0_channel', 'adc1_channel'):
        message.append(int(arduino_config[key]))
    send_to_arduino(message)
    # print("Gain in ADC initialisation is:", receive_from_arduino())
    if receive_from_arduino() == PC_OK:
        print("ADC Initalisation: DONE!")
        print("-----------------------")
    else:
        identify_error()
        arduino_port.close()
        raise IOError("ADC Initalisation: FAILED!")

# ======================================


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

    # if sys.platform.startswith('win'):
        # # ports = list(serial.tools.list_ports.comports())
        # ports = qts.QSerialPortInfo.availablePorts()
        # for port in ports:
            # if 'Arduino Micro' in port[1]: results.append(port[0])
    # elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # # this excludes your current terminal "/dev/tty"
        # ports = glob.glob('/dev/tty[A-Za-z]*')
    # elif sys.platform.startswith('darwin'):
        # ports = glob.glob('/dev/tty[A-Za-z]*')
    # else:
        # raise EnvironmentError('Unsupported platform')

    arduino_ports = [port
                     for port in qts.QSerialPortInfo.availablePorts()
                     if 'Arduino Micro' in port.description()]
    if not arduino_ports:
        raise IOError("Couldn't find port to Arduino!")

    for port_info in arduino_ports:
        port = qts.QSerialPort(port_info)
        if not port.open(port.ReadWrite):
            continue
        port.close()
        return port_info

    raise IOError("Unable to connect to Arduino, try "
                  "disconnecting and connecting again.")

# ======================================
def set_voltage_and_measure(energy, settle_time = 0):
    """
    Function initialise dac and sends dac value to arduino.

    Parameter
    ----------
    energy: float
         Nominal electron energy in electronvolts to be set
    settle_time: int
    """
    v_ref_dac = config.getfloat('measurement_settings', 'v_ref_dac')
    PC_SET_VOLTAGE = config.getint('communication_bytes', 'PC_SET_VOLTAGE')
    dac_out_vs_nominal_energy = 10/1000  # 10V / 1000 eV
    ouput_gain = 4  # Gain of the output stage on board
    conversion_factor = dac_out_vs_nominal_energy * 65536 / (v_ref_dac * ouput_gain)

    # TODO: nominal energy to actual energy after calibration
    dac_value_temp = int(round(energy * conversion_factor))
    if dac_value_temp >= 65536:
        dac_value_temp = 65535
    if dac_value_temp <= 0:
        dac_value_temp = 0
    dac_send = dac_value_temp.to_bytes(2, sign) + settle_time.to_bytes(2, sign)
    # print("DAC_send = " ,dac_value_temp)
    # print("DAC_tobytes =" ,dac_value_temp.to_bytes(2, sign), " t= ",round(1000*time.time()))
    send_to_arduino(PC_SET_VOLTAGE)
    send_to_arduino(dac_send)
    # dac_send = send_to_arduino(dac_send)
    if receive_from_arduino() == PC_OK:
        print("Voltage is set, triggering")
    else:
        identify_error()
        arduino_port.close()
# ======================================

def start_autogain():
    """
    Function starts autogain-function on arduino. Autogain-funtion measures
    adc value and sets gain, that measured value is bigger than 25% of the
    measurment range. Starts with gain 1 (no gain). Maximum gain is 128
    """
    PC_AUTOGAIN = config.getint('communication_bytes', 'PC_AUTOGAIN')
    send_to_arduino(PC_AUTOGAIN)
    print("Autogain: IN PROCESS...")
    if receive_from_arduino() == PC_OK:
        # print("Gain found =",pow(2,int(receive_from_arduino())))#for debugging, can be deleted
        # print("Gain found =", int(receive_from_arduino()))
        pass
    else:
        identify_error()
        arduino_port.close()
        raise IOError("Arduino doesnt react when python tries to initialise "
                      "the Autogain Function, Program stops")

# ======================================

def create_csv(multilist, path):
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
    # timestamp = datetime.now()
    df = pd.DataFrame(multilist, columns=['ADC_VAL', 'DAC_Val', 'DeltaTime'])
    print("-----------------------")
    print(df)
    print("-----------------------")
    df.to_csv(os.path.join(path, "sample%s.csv" % i), index=True)
    print("CSV File was created")

# ======================================

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
        max_energy = 4*config.getfloat('measurement_settings', 'v_ref_dac')
        start_energy = -1
        while start_energy < 0 or start_energy > max_energy:
            start_energy = float(input(
                "\nPlease enter the energy to START with in Volt. "
                f"Allowed are values between 0 V and "
                f"{max_energy*100:.1f} V:\n"
                ))/100
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
            delta_energy = float(
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
    while os.path.exists(os.path.join(zip_path, str(i), ".zip")):
        i += 1

    with ZipFile(os.path.join(zip_path, str(i), ".zip")) as zip_file:
        for f in fnames:
            zip_file.write(f, os.path.basename(f))

def calibrate_adcs(arduino_config):
    """
    Send PC_CALIBRATION to arduino which does the calibration for all
    gains for the selected channels and saves the values for later use.

    Parameter
    ----------
    arduino_config: dict
    """
    PC_CALIBRATION = config.getint('communication_bytes', 'PC_CALIBRATION')
    send_to_arduino(PC_CALIBRATION)
    message = []
    for key in ('update_rate', 'adc0_channel', 'adc1_channel'):
        message.append(int(arduino_config[key]))
    send_to_arduino(message)
    if receive_from_arduino() == PC_OK:
        print("ADC calibration: DONE!")
        print("-----------------------")
    else:
        identify_error()
        arduino_port.close()
        raise IOError("ADC calibration: FAILED!")


def identify_error():
    """
    Is called upon when PC_ERROR is returned from the Arduino to the
    PC. It compares the trace back byte to the ones saved in the config
    and prints the according key
    """
    arduino_states = config['arduino_states']
    error_bytes = config['error_bytes']
    ErrorMessage = receive_from_arduino()
    for key, value in arduino_states.items():
        if int(value) == ErrorMessage[0]:
            print(key)
    for key, value in error_bytes.items():
        if int(value) == ErrorMessage[1]:
            print(key)

def prepare_for_measurement(configuration_section):

    adc_config = config[configuration_section]
    dac_first_settletime = config.getint('measurement_settings',
                                         'dac_first_settle_time')
    start_energy = config.getfloat('measurement_settings', 'start_energy')

    dac_energy = energy_calibration_curve()

    # Connect with Arduino, reset it and flush input buffer
    portname = serial_ports()
    arduino_port = connect_to_arduino(portname)
    arduino_port.flushInput()

    # Calibrate and initialize ADCs, DAC and auto-gain,
    # then re-initialize the ADC with the right gain
    calibrate_adcs(adc_config)
    initialise_adcs(adc_config)
    set_voltage_and_measure(dac_energy(start_energy), dac_first_settletime)
    receive_from_arduino()
    receive_from_arduino()
    receive_from_arduino()
    start_autogain()

def measure_iv_video():
    prepare_for_measurement('iv_movie_configuration')
    PC_RESET = config.getint('communication_bytes', 'PC_RESET')
    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    path = config['measurement_settings']['path']
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    camera_settle_time = config.getint('camera_settings', 'camera_settle_time')
    live_mode = (config['camera_settings']['live_mode'] == 'True')
    I0_conversion_factor = config.getint('leed_hardware',
                                         'I0_conversion_factor')

    start_energy, end_energy, delta_energy = TUI()
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    delta_energy = config.getfloat('measurement_settings', 'delta_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')

    # Connect with Camera and initialize it
    CameraObject = Camera.camera_from_ini(
        config['camera_settings']['class_name'])
    CameraObject.initialize(config['camera_settings'])

    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()

    npoints = int((end_energy-start_energy)/delta_energy)
    if npoints <= 0:
        print("Please check the starting and ending energy and the stepsize.")
        raise RuntimeError("Number of steps could not be calculated!")
    nominal_energy_csv = np.linspace(start_energy, end_energy, npoints+1)

    # This keyboard interrupt thing will not exist. Process interruption will
    # be handled with a signal
    # TODO: I have no clue what this does but I removed the call for a measurement as set_voltage_and_measure already does that. We need to redo this.
    try:
        for actual_energy in nominal_energy_csv:
            print("Energy:%.2f V" % actual_energy)
            # if camera is used in "exposure" mode
            while not live_mode:
                if not CameraObject.callback_data.dac_busy:
                    set_voltage_and_measure(actual_energy, dac_settle_time)
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
            # print('Gain: %i \nDAC-Value: %f' % (int(receive_from_arduino()),
            #                                   actual_energy))
            print('nDAC-Value: %f' % actual_energy)
            adc0_value = bytes_to_float(receive_from_arduino())
            adc1_value = bytes_to_float(receive_from_arduino())
            adc2_value = bytes_to_float(receive_from_arduino())
            # adc0 is ADC0, adc1 is ADC1 and adc2 is the LM35 in the Arduino code
            print("ADC0_VAlUE:", adc0_value*I0_conversion_factor)
            print("ADC1_VAlUE:", adc1_value)
            print("ADC2_VAlUE:", adc2_value,
                  " t = ", round(1000*time.time()))
            print("#########################")
            adc0_value_csv.append(adc0_value*I0_conversion_factor)
            adc1_value_csv.append(adc1_value)
            adc2_value_csv.append(adc2_value)
            time_temp = time.time() - timestart
            timestamp.append(time_temp)

    # finally:
    except KeyboardInterrupt:
        CameraObject.stop_camera()
        # send_to_arduino(PC_RESET)
        # Reset would clear calibration: not necessary
        set_voltage_and_measure(dac_value_end)
        arduino_port.close()

    time.sleep(3)  # To process the last frame, the callback needs some time
    CameraObject.stop_camera()
    set_voltage_and_measure(dac_value_end)
    create_csv(list(zip(adc0_value_csv, adc1_value_csv, adc2_value_csv,
                        nominal_energy_csv, timestamp)), path)
    # before terminating, pack all the necessary data into a single zip file
    # pack_measurements()

def test_iv_video():
    prepare_for_measurement('iv_movie_configuration')
    # PC_RESET = config.getint('communication_bytes', 'PC_RESET')
    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    I0_conversion_factor = config.getint('leed_hardware',
                                         'I0_conversion_factor')

    start_energy = config.getfloat('measurement_settings', 'start_energy')
    delta_energy = config.getfloat('measurement_settings', 'delta_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')

    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()

    npoints = int((end_energy-start_energy)/delta_energy)
    if npoints <= 0:
        print("Please check the starting and ending energy and the stepsize.")
        raise RuntimeError("Number of steps could not be calculated!")
    nominal_energy_csv = np.linspace(start_energy, end_energy, npoints+1)

    for actual_energy in nominal_energy_csv:
        set_voltage_and_measure(actual_energy, dac_settle_time)
        print(actual_energy)
        data = receive_from_arduino()
        if len(data) == 1:
            identify_error()
            break
        else:
            adc0_value = bytes_to_float(data)
        adc1_value = bytes_to_float(receive_from_arduino())
        adc2_value = bytes_to_float(receive_from_arduino())
        print("ADC0_VAlUE:", adc0_value*I0_conversion_factor)
        print("ADC1_VAlUE:", adc1_value)
        print("ADC2_VAlUE:", adc2_value,
              " t = ", round(1000*time.time()))
        print("#########################")
        adc0_value_csv.append(adc0_value*I0_conversion_factor)
        adc1_value_csv.append(adc1_value)
        adc2_value_csv.append(adc2_value)
        time_temp = time.time() - timestart
        timestamp.append(time_temp)

    set_voltage_and_measure(dac_value_end)
    fit_polynomial = Polynomial.fit(nominal_energy_csv, adc0_value_csv, 5)
    print(fit_polynomial, flush=True)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_ylabel('adc0_value')
    ax1.set_xlabel('energy')
    ax1.plot(nominal_energy_csv, adc0_value_csv, '.')
    ax1.plot(nominal_energy_csv, fit_polynomial(nominal_energy_csv), '-')

    ax2.set_ylabel('residuals')
    ax2.set_xlabel('true energy')
    ax2.plot(nominal_energy_csv,
             np.subtract(fit_polynomial(nominal_energy_csv), adc0_value_csv),
             '.')
    plt.show()

def change_measurement_mode(measurement_mode):
    """
    Send measurement mode to arduino, which decides if the arduino will do
    continuous measurememts or not.
    """
    PC_CHANGE_MEASUREMENT_MODE = config.getint('communication_bytes',
                                               'PC_CHANGE_MEAS_MODE')
    PC_ERROR = config.getint('communication_bytes',
                             'PC_ERROR').to_bytes(1, sign)
    send_to_arduino(PC_CHANGE_MEASUREMENT_MODE)
    message = []
    use_this_mode = config.getint('measurement_settings',
                                  measurement_mode)
    message.append(use_this_mode)
    measurement_time = config.getint('measurement_settings',
                                     'measurement_time').to_bytes(2, sign)
    message.append(measurement_time[0])
    message.append(measurement_time[1])
    # print(message)
    send_to_arduino(message)
    msg = receive_from_arduino()
    while (msg != PC_ERROR) and (msg != PC_OK):
        msg = receive_from_arduino()
    if msg == PC_OK:
        print("Measurement mode has been set.")
    else:
        identify_error()
        arduino_port.close()
        raise IOError("Setting the measurement mode failed.")

def quick_measurements(measure_what='I0', settle_time_scaling=1):
    """"""
    conversion_factor = 1
    if measure_what == 'I0':
        config_section = 'iv_movie_configuration'
        conversion_factor = config.getint('leed_hardware',
                                          'I0_conversion_factor')
        y_title = 'I0 (uA)'
    elif measure_what == 'HV':
        config_section = 'measure_filament_configuration'
        y_title = 'Energy (eV)'
    else:
        raise ValueError("Unknown measurement mode")
    dac_energy = energy_calibration_curve()
    prepare_for_measurement(config_section)
    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    dac_first_settle_time = config.getint('measurement_settings',
                                          'dac_first_settle_time')
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    continuous_measurement_points = config.getint('measurement_settings',
                                                 'continuous_measurement_points')
    continuous_measurement_loops = config.getint('measurement_settings',
                                                'continuous_measurement_loops')
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')
    update_rate_raw = config.getint(config_section, 'update_rate')
    update_rate = config.getint('adc_update_rate', str(update_rate_raw))

    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()
    for _ in range(3):
        set_voltage_and_measure(end_energy, dac_settle_time)
        receive_from_arduino()
        receive_from_arduino()
        receive_from_arduino()
        set_voltage_and_measure(start_energy, dac_first_settle_time)
    j = 1
    while j <= continuous_measurement_loops:
        set_voltage_and_measure(start_energy, dac_first_settle_time)
        change_measurement_mode('continuous_measurement_yes')
        print(dac_settle_time)
        set_voltage_and_measure(end_energy, dac_settle_time)

        i = 1
        adc0_value_this_loop = []
        while i <= continuous_measurement_points:
            i += 1
            adc0_value = bytes_to_float(receive_from_arduino())
            adc1_value = bytes_to_float(receive_from_arduino())
            adc2_value = bytes_to_float(receive_from_arduino())
            # print("ADC0_VAlUE:", adc0_value*conversion_factor)
            # print("ADC1_VAlUE:", adc1_value)
            # print("ADC2_VAlUE:", adc2_value,
                  # " t = ", round(1000*time.time()))
            # print("#########################")
            adc0_value_this_loop.append(adc0_value*conversion_factor)
            adc1_value_csv.append(adc1_value)
            adc2_value_csv.append(adc2_value)
            time_temp = time.time() - timestart
            timestamp.append(time_temp)
        dac_settle_time *= settle_time_scaling
        adc0_value_csv.append(adc0_value_this_loop)
        j += 1

        change_measurement_mode('continuous_measurement_no')

    set_voltage_and_measure(dac_value_end)

    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    times = np.arange(continuous_measurement_points)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_ylabel(y_title)
    ax1.set_xlabel('Time since set voltage (ms)')
    averaged = np.zeros(np.shape(times))
    for i, data in enumerate(adc0_value_csv):
        averaged = np.add(averaged, data)
        ax1.plot((times + 3)*1000/update_rate + dac_settle_time*settle_time_scaling**i,
                 data, '.')

    ax2.set_ylabel('Average ' + y_title)
    ax2.set_xlabel('Time since set voltage (ms)')
    ax2.plot((times + 3)*1000/update_rate,
             averaged/continuous_measurement_loops,
             '.')
    plt.show()
    return adc0_value_csv

def test_single_ramp_quick_measurements(measure_what='I0'):
    conversion_factor = 1
    if measure_what == 'I0':
        config_section = 'iv_movie_configuration'
        conversion_factor = config.getint('leed_hardware',
                                          'I0_conversion_factor')
        y_title = 'I0 (uA)'
    elif measure_what == 'HV':
        config_section = 'measure_filament_configuration'
        y_title = 'Energy (eV)'
    else:
        raise ValueError("Unknown measurement mode")

    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    dac_first_settle_time = config.getint('measurement_settings',
                                          'dac_first_settle_time')
    continuous_measurement_points = config.getint('measurement_settings',
                                                 'continuous_measurement_points')
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    delta_energy = config.getfloat('measurement_settings', 'delta_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')

    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()

    npoints = int((end_energy-start_energy)/delta_energy)
    if npoints <= 0:
        print("Please check the starting and ending energy and the stepsize.")
        raise RuntimeError("Number of steps could not be calculated!")
    nominal_energy_csv = np.linspace(start_energy, end_energy, npoints+1)

    for actual_energy in nominal_energy_csv:
        change_measurement_mode('continuous_measurement_yes')
        set_voltage_and_measure(actual_energy, 0)
        i = 1
        adc0_value_this_loop = []
        while i <= continuous_measurement_points:
            i += 1
            adc0_value = bytes_to_float(receive_from_arduino())
            adc1_value = bytes_to_float(receive_from_arduino())
            adc2_value = bytes_to_float(receive_from_arduino())
            adc0_value_this_loop.append(adc0_value*conversion_factor)
            adc1_value_csv.append(adc1_value)
            adc2_value_csv.append(adc2_value)
            time_temp = time.time() - timestart
            timestamp.append(time_temp)
        adc0_value_csv.append(adc0_value_this_loop)

        change_measurement_mode('continuous_measurement_no')

    set_voltage_and_measure(start_energy, dac_first_settle_time)

    return adc0_value_csv

def multiple_ramps_quick_measurements(measure_what='I0'):
    conversion_factor = 1
    if measure_what == 'I0':
        config_section = 'iv_movie_configuration'
        conversion_factor = config.getint('leed_hardware',
                                          'I0_conversion_factor')
        y_title = 'I0 (uA)'
    elif measure_what == 'HV':
        config_section = 'measure_filament_configuration'
        y_title = 'Energy (eV)'
    else:
        raise ValueError("Unknown measurement mode")
    prepare_for_measurement(config_section)
    update_rate_raw = config.getint(config_section, 'update_rate')
    update_rate = config.getint('adc_update_rate', str(update_rate_raw))
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    num_ramps = config.getint('measurement_settings',
                              'continuous_measurement_loops')
    continuous_measurement_points = config.getint('measurement_settings',
                                                 'continuous_measurement_points')
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')
    adc0_value_csv = []
    i=1
    while i <= num_ramps:
        adc0_value_csv.append(test_single_ramp_quick_measurements(measure_what))
        i += 1
    # times in milliseconds. The +3 is because of triggering.
    n_points_ave = round(0.25*continuous_measurement_points)
    times = (np.arange(continuous_measurement_points) + 3)*1000/update_rate
    adc0_values = np.asarray(adc0_value_csv)
    average = np.sum(adc0_values, axis=0)/num_ramps
    shift = np.sum(average[:,-n_points_ave:], axis=1)/n_points_ave
    steps = (average.T - shift).T
    standard_deviation = np.std(steps[:,-n_points_ave:], axis=1, ddof=1)
    standard_dev_min = min(dev for dev in standard_deviation if dev > 0)
    num_voltage_steps = len(standard_deviation)

    for j, stdev in enumerate(standard_deviation):
        if (stdev == 0):
            standard_deviation[j] = standard_dev_min
        print(standard_deviation[j])

    relaxation_time = []
    sigma_mult = 2.8
    for step, stdev in zip(steps, standard_deviation):
        i = 1
        while i < continuous_measurement_points - 2:
            if (abs(step[-i]) > sigma_mult*stdev
                and abs(step[-i-1]) > sigma_mult*stdev
                    and abs(step[-i-2]) > sigma_mult*stdev):
                break
            i += 1
        relaxation_time.append(continuous_measurement_points - i)

    relaxation_time = (np.asarray(relaxation_time) + 3)*1000/update_rate

    energy_scale = np.linspace(start_energy, end_energy, num_voltage_steps)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    ax1.set_ylabel(y_title)
    ax1.set_xlabel('Time since set voltage (ms)')
    for sigma, step in zip(standard_deviation, steps):
        _ = ax1.plot(times, step, '.')[0]
        ax1.plot((times[0], times[-1]), (sigma_mult*sigma, sigma_mult*sigma),
                 '-', color=_.get_color())
        ax1.plot((times[0], times[-1]), (-sigma_mult*sigma, -sigma_mult*sigma),
                 '-', color=_.get_color())
    ax2.set_ylabel("Relaxation time (ms)")
    ax2.set_xlabel("Energy (eV)")
    ax2.plot(energy_scale, relaxation_time, '.')
    ax3.set_ylabel("Sigma (eV)")
    ax3.set_xlabel("Energy (eV)")
    ax3.plot(energy_scale, standard_deviation, '.')
    plt.show()
    set_voltage_and_measure(dac_value_end)
    return times, steps

def long_term_measurements(measure_what='I0', interval=5):
    conversion_factor = 1
    if measure_what == 'I0':
        config_section = 'iv_movie_configuration'
        conversion_factor = config.getint('leed_hardware',
                                          'I0_conversion_factor')
        y_title = 'I0 (uA)'
    elif measure_what == 'HV':
        config_section = 'measure_filament_configuration'
        y_title = 'Energy (eV)'
    else:
        raise ValueError("Unknown measurement mode")
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    prepare_for_measurement(config_section)

    counter = 720
    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()

    while counter:
        adc_measure_only()
        adc0_value = bytes_to_float(receive_from_arduino())
        adc1_value = bytes_to_float(receive_from_arduino())
        adc2_value = bytes_to_float(receive_from_arduino())
        adc0_value_csv.append(adc0_value*conversion_factor)
        adc1_value_csv.append(adc1_value)
        adc2_value_csv.append(adc2_value)
        time_temp = time.time() - timestart
        timestamp.append(time_temp)
        counter -= 1
        time.sleep(interval)

    set_voltage_and_measure(dac_value_end)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    ax1.set_ylabel(y_title)
    ax1.set_xlabel('Time since set voltage (s)')
    ax1.plot(timestamp, adc0_value_csv, '.')
    ax2.set_ylabel("ADC1 (V)")
    ax2.set_xlabel("Time since set voltage (s)")
    ax2.plot(timestamp, adc1_value_csv, '.')
    ax3.set_ylabel("Temperature (°C)")
    ax3.set_xlabel("Time since set voltage (s)")
    ax3.plot(timestamp, adc2_value_csv, '.')
    plt.show()
    return timestamp, adc0_value_csv, adc1_value_csv, adc2_value_csv

def quick_up_down_measurements(measure_what='I0', settle_time_scaling=1):
    """"""
    conversion_factor = 1
    if measure_what == 'I0':
        config_section = 'iv_movie_configuration'
        conversion_factor = config.getint('leed_hardware',
                                          'I0_conversion_factor')
        y_title = 'I0 (uA)'
    elif measure_what == 'HV':
        config_section = 'measure_filament_configuration'
        y_title = 'Energy (eV)'
    else:
        raise ValueError("Unknown measurement mode")
    dac_energy = energy_calibration_curve()
    prepare_for_measurement(config_section)
    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    dac_first_settle_time = config.getint('measurement_settings',
                                         'dac_first_settle_time')
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    continuous_measurement_points = config.getint('measurement_settings',
                                                 'continuous_measurement_points')
    continuous_measurement_loops = config.getint('measurement_settings',
                                                 'continuous_measurement_loops')
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')
    update_rate_raw = config.getint(config_section, 'update_rate')
    update_rate = config.getint('adc_update_rate', str(update_rate_raw))

    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()
    for _ in range(3):
        set_voltage_and_measure(end_energy, dac_settle_time)
        receive_from_arduino()
        receive_from_arduino()
        receive_from_arduino()
        set_voltage_and_measure(start_energy, dac_first_settle_time)
    j = 1
    while j <= continuous_measurement_loops:
        set_voltage_and_measure(start_energy, dac_first_settle_time)
        change_measurement_mode('continuous_measurement_yes')
        set_voltage_and_measure(end_energy+0.1*np.sign(end_energy-start_energy), 0)
        set_voltage_and_measure(end_energy, dac_settle_time)

        i = 1
        adc0_value_this_loop = []
        while i <= continuous_measurement_points:
            i += 1
            adc0_value = bytes_to_float(receive_from_arduino())
            adc1_value = bytes_to_float(receive_from_arduino())
            adc2_value = bytes_to_float(receive_from_arduino())
            adc0_value_this_loop.append(adc0_value*conversion_factor)
            adc1_value_csv.append(adc1_value)
            adc2_value_csv.append(adc2_value)
            time_temp = time.time() - timestart
            timestamp.append(time_temp)
        dac_settle_time *= settle_time_scaling
        adc0_value_csv.append(adc0_value_this_loop)
        j += 1

        change_measurement_mode('continuous_measurement_no')

    set_voltage_and_measure(dac_value_end)

    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    times = np.arange(continuous_measurement_points)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_ylabel(y_title)
    ax1.set_xlabel('Time since set voltage (ms)')
    averaged = np.zeros(np.shape(times))
    for i, data in enumerate(adc0_value_csv):
        averaged = np.add(averaged, data)
        ax1.plot((times + 3)*1000/update_rate + dac_settle_time*settle_time_scaling**i,
                 data, '.')

    ax2.set_ylabel('Average ' + y_title)
    ax2.set_xlabel('Time since set voltage (ms)')
    ax2.plot((times + 3)*1000/update_rate,
             averaged/continuous_measurement_loops,
             '.')
    plt.show()
    plt.yscale('log')
    yf = rfft(averaged[:]/continuous_measurement_loops)
    xf = rfftfreq(continuous_measurement_points, 1 / update_rate)
    plt.plot(xf, np.abs(yf))
    plt.show()

    return adc0_value_csv

def quick_up_down_measurement_ramp(measure_what='I0', settle_time_scaling=1):
    """"""
    conversion_factor = 1
    if measure_what == 'I0':
        conversion_factor = config.getint('leed_hardware',
                                          'I0_conversion_factor')
    elif measure_what == 'HV':
        pass
    else:
        raise ValueError("Unknown measurement mode")
    dac_energy = energy_calibration_curve()
    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    dac_first_settle_time = config.getint('measurement_settings',
                                         'dac_first_settle_time')
    continuous_measurement_points = config.getint('measurement_settings',
                                                 'continuous_measurement_points')
    continuous_measurement_loops = config.getint('measurement_settings',
                                                 'continuous_measurement_loops')
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')
    delta_energy = config.getfloat('measurement_settings', 'delta_energy')

    adc0_value_csv = []
    adc1_value_csv = []
    adc2_value_csv = []
    timestamp = []
    timestart = time.time()

    npoints = int((end_energy-start_energy)/delta_energy)
    if npoints <= 0:
        print("Please check the starting and ending energy and the stepsize.")
        raise RuntimeError("Number of steps could not be calculated!")
    nominal_energy_csv = np.linspace(start_energy, end_energy, npoints+1)

    if start_energy > 0:
        set_voltage_and_measure(start_energy - delta_energy,
                                dac_first_settle_time)
    else:
        set_voltage_and_measure(0, dac_first_settle_time)

    for actual_energy in nominal_energy_csv:
        adc0_value_this_loop = []
        set_voltage_and_measure(actual_energy+0.1*np.sign(delta_energy), 0)
        change_measurement_mode('continuous_measurement_yes')
        set_voltage_and_measure(actual_energy, dac_settle_time)
        i = 1
        while i <= continuous_measurement_points:
            i += 1
            adc0_value = bytes_to_float(receive_from_arduino())
            adc1_value = bytes_to_float(receive_from_arduino())
            adc2_value = bytes_to_float(receive_from_arduino())
            adc0_value_this_loop.append(adc0_value*conversion_factor)
            adc1_value_csv.append(adc1_value)
            adc2_value_csv.append(adc2_value)
            time_temp = time.time() - timestart
            timestamp.append(time_temp)
        dac_settle_time *= settle_time_scaling
        adc0_value_csv.append(adc0_value_this_loop)
        change_measurement_mode('continuous_measurement_no')

    return adc0_value_csv

def multiple_up_down_measurement_ramps(measure_what='I0'):
    if measure_what == 'I0':
        config_section = 'iv_movie_configuration'
        y_title = 'I0 (uA)'
    elif measure_what == 'HV':
        config_section = 'measure_filament_configuration'
        y_title = 'Energy (eV)'
    else:
        raise ValueError("Unknown measurement mode")
    prepare_for_measurement(config_section)
    update_rate_raw = config.getint(config_section, 'update_rate')
    update_rate = config.getint('adc_update_rate', str(update_rate_raw))
    num_ramps = config.getint('measurement_settings',
                              'continuous_measurement_loops')
    continuous_measurement_points = config.getint('measurement_settings',
                                                 'continuous_measurement_points')
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    dac_first_settle_time = config.getint('measurement_settings',
                                         'dac_first_settle_time')
    adc0_value_csv = []
    i=1
    while i <= num_ramps:
        adc0_value_csv.append(quick_up_down_measurement_ramp(measure_what))
        i += 1
    set_voltage_and_measure(dac_value_end, dac_first_settle_time)
    # times in milliseconds. The +3 is because of triggering.
    n_points_ave = round(0.25*continuous_measurement_points)
    times = (np.arange(continuous_measurement_points) + 3)*1000/update_rate
    adc0_values = np.asarray(adc0_value_csv)
    average = np.sum(adc0_values, axis=0)/num_ramps
    shift = np.sum(average[:,-n_points_ave:], axis=1)/n_points_ave
    steps = (average.T - shift).T

    average_heights = []
    average_heights = shift[1:] - shift[:-1]
    average_heights = np.insert(average_heights, 0, shift[0])

    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_ylabel(y_title)
    ax1.set_xlabel('Time since set voltage (ms)')
    xf = rfftfreq(continuous_measurement_points, 1 / update_rate)
    ax2.set_yscale('log')
    ax2.set_ylabel('FFT (absolute value)')
    ax2.set_xlabel('Frequency (Hz)')
    for step in steps:
        ax1.plot(times, step, '.')
        ax1.plot((times[0], times[-1]), (average_heights, average_heights), '-')
        yf = rfft(step/num_ramps)
        ax2.plot(xf, np.abs(yf))
    plt.show()
    return times, steps

def main():
    global arduino_port
    app = qtw.QApplication(sys.argv)
    # measure_iv_video()
    # test_iv_video()
    # ret_value = quick_measurements('I0', 2)
    # ret_value = quick_measurements('HV', 2)
    # calibrate_real_energy_scale()
    # ret_value = multiple_ramps_quick_measurements('HV')
    # ret_value = long_term_measurements('HV')
    # ret_value = long_term_measurements('I0')
    ret_value = quick_up_down_measurements('I0', 2)
    # ret_value = multiple_up_down_measurement_ramps('I0')
    arduino_port.close()
    sys.exit(app.exec_())
    print('Arduino Disconnected')
    return ret_value


def do_nothing_on_purpose():  # MR: just a mock I needed for my presentation

    start_energy, end_energy, delta_energy = TUI()
    actual_energy = start_energy

    time.sleep(1) # time to settle connection
    print("Arduino Initialisation: DONE!")

    time.sleep(0.5) # time to settle connection
    print("ADC Initalisation: DONE!")
    print("-----------------------")
    time.sleep(0.2) # time to settle connection
    print("Autogain: IN PROCESS...")

    time.sleep(0.2) # time to settle connection
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

def calibrate_real_energy_scale():
    """
    This function measures the voltage on the filament that we get in
    return for the requested voltage. After this the collected data is
    used to do a polynomial fit which can be used to obtain the correct
    voltage on the filament.
    """
    prepare_for_measurement('measure_filament_configuration')
    dac_settle_time = config.getint('measurement_settings', 'dac_settle_time')
    # start_energy, end_energy, delta_energy = TUI()
    start_energy = config.getfloat('measurement_settings', 'start_energy')
    delta_energy = config.getfloat('measurement_settings', 'delta_energy')
    end_energy = config.getfloat('measurement_settings', 'end_energy')
    dac_value_end = config.getfloat('measurement_settings', 'dac_value_end')
    true_energy_csv = []

    npoints = int((end_energy-start_energy)/delta_energy)
    if npoints <= 0:
        print("Please check the starting and ending energy and the stepsize.")
        raise RuntimeError("Number of steps could not be calculated!")
    nominal_energy_csv = np.linspace(start_energy, end_energy, npoints+1)

    for actual_energy in nominal_energy_csv:
        set_voltage_and_measure(actual_energy, dac_settle_time)
        data = receive_from_arduino()
        if len(data) == 1:
           identify_error()
           break
        else:
            adc0_value = bytes_to_float(data)
        receive_from_arduino()
        receive_from_arduino()
        print("ADC0_VAlUE:", adc0_value)
        print("#########################")
        true_energy_csv.append(adc0_value)

    true_energy_csv = np.asarray(true_energy_csv)

    set_voltage_and_measure(dac_value_end)
    fit_polynomial = Polynomial.fit(true_energy_csv, nominal_energy_csv, 1, domain=[-10,1100], window=[-10, 1100])

    print(fit_polynomial.domain)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_ylabel('fit_polynomial')
    ax1.set_xlabel('true energy')
    ax2.set_ylabel('residuals')
    ax2.set_xlabel('true energy')
    print("[", ", ".join(f"{ci:.20f}" for ci in fit_polynomial.coef), "]", sep="")

    ax1.plot(true_energy_csv, nominal_energy_csv, '.')
    ax1.plot(true_energy_csv, fit_polynomial(true_energy_csv), '-')
    ax2.plot(true_energy_csv,
             np.subtract(fit_polynomial(true_energy_csv), nominal_energy_csv),
             '.')
    #plt.ylim(0, 10)
    plt.show()

def energy_calibration_curve():
    calibration_coef = ast.literal_eval(
        config['energy_calibration']['coefficients']
        )
    calibration_domain = ast.literal_eval(
        config['energy_calibration']['domain']
        )
    dac_energy = Polynomial(calibration_coef, domain=calibration_domain,
                            window=calibration_domain)
    return dac_energy

if __name__ == '__main__':
    t1 = time.time()
    for i in range(1):
        main()
        print('#########################')
        print('DONE!')
        print("Time for whole procedure:", (time.time() - t1))

