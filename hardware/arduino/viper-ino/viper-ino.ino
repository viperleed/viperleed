/* ViPErLEED - Firmware for Arduino hardware controller

@author: Bernhard Mayr
@author: Michael Schmid (@schmid-iap)
@author: Michele Riva (@michele-riva)
@author: Florian Dörr (@FlorianDoerr)
Date: 09.02.2022
*/

#include "viper-ino.h"  // Arduino-related constants, globals and functions.
                        // Includes serial, states, ADC and DAC.

#define DEBUG   false   // Debug mode, writes to serial line,
                        // for use in serial monitor.

// The box ID is an indentifier that is necessary for the PC to know what
// type of Arduino it is handling. 1 is the identifier of a ViPErino
// controller that performs LEED measurements. The box ID must never be
// changed!
#define BOX_ID  1

// Firmware version (MAX: v255.255). CURENTLY: v0.10
#define FIRMWARE_VERSION_MAJOR    0  // max 255
#define FIRMWARE_VERSION_MINOR    10  // max 255


/** ---------------------------- INITIALIZATION ---------------------------- **/

void setup() {
    /** Set up the Arduino board ans all the hardware.

    The setup routine runs once on power-up or on hardware reset.
    */
    #if DEBUG
        // initialize serial communication (emulation
        // on USB) at 9600 bits per second for debug
        delay(1000);
        Serial.begin(9600);
        delay(1000);
    #endif

    analogReference(INTERNAL); // Used for LM35 and for checking relay presence

    // Set all inverted chip-select pins to high via pull-up
    // (to avoid low glitch)
    setChipSelectHigh(CS_DAC);
    setChipSelectHigh(CS_ADC_0);
    setChipSelectHigh(CS_ADC_1);
    pinMode(SCK, OUTPUT);  // Should not be needed, but it did not work without
                           // Probably just because SPI.begin call is missing?
    pinMode(MOSI, OUTPUT);

    // Reset ADC and DAC, set DAC output to zero,
    // and set variables to well-defined defaults
    reset();
}


/** ----------------------------- STATE MACHINE ---------------------------- **/

void loop() {
    /**Read from serial line and handle the state machine.**/
    // See if there is anything coming from the PC.
    readFromSerial();

    // Then see if the current state needs updating,
    // depending on the message received.
    updateState();

    // Finally do what's needed in the currentState.
    switch (currentState){
        case STATE_GET_CONFIGURATION:
            getConfiguration();
            break;
        case STATE_CALIBRATE_ADCS:
            calibrateADCsAtAllGains();
            break;
        case STATE_SET_UP_ADCS:
            prepareADCsForMeasurement();
            break;
        case STATE_SET_VOLTAGE:
            setVoltageWaitAndTrigger();
            break;
        case STATE_MEASURE_ADCS:
            measureADCs();
            break;
        case STATE_ADC_VALUES_READY:
            sendMeasuredValues();
            break;
        case STATE_AUTOGAIN_ADCS:
            findOptimalADCGains();
            break;
        case STATE_ERROR:
            handleErrors();
            break;
        case STATE_CHANGE_MEASUREMENT_MODE:
            changeMeasurementMode();
            break;
        case STATE_SET_SERIAL_NR:
            setSerialNr();
            break;
    }
}


void updateState() {
    /**
    If there is an unprocessed serial message, decide whether
    this requires us to change the state of the Arduino, and
    do some preparation that may be needed before entering
    the new state. All data messages have to be at least two
    bytes long. Messages that are one byte long are assumed
    to be commands.

    Reads
    -----
    newMessage, data_received

    Writes
    ------
    initialTime, currentState, waitingForDataFromPC, calibrationGain

    Goes to state
    -------------
    May go anywhere
    **/
    if (not newMessage) return;

    if (msgLength > 1) return; // We received data, not a command
    // Note that all data messages have to be at least two bytes long.

    if (data_received[0] != PC_CONFIGURATION
        and data_received[0] != PC_RESET
        and data_received[0] != PC_STOP
        and data_received[0] != PC_CHANGE_MEAS_MODE
        and data_received[0] != PC_SET_VOLTAGE_ONLY
        and data_received[0] != PC_SET_SERIAL_NR
        and hardwareNotKnown()) return;

    switch(data_received[0]){
        case PC_CONFIGURATION:
            initialTime = millis();
            currentState = STATE_GET_CONFIGURATION;
            break;
        case PC_CALIBRATION:
            waitingForDataFromPC = true;
            initialTime = millis();
            calibrationGain = 0;
            currentState = STATE_CALIBRATE_ADCS;
            break;
        case PC_SET_UP_ADCS:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_UP_ADCS;
            break;
        case PC_SET_VOLTAGE:
            takeMeasurements = true;
            waitingForDataFromPC = true;
            initialTime = millis();
            nextVoltageStep = 0;
            currentState = STATE_SET_VOLTAGE;
            break;
        case PC_AUTOGAIN:
            encodeAndSend(PC_OK);
            initialTime = millis();
            prepareForAutogain();
            currentState = STATE_AUTOGAIN_ADCS;
            break;
        case PC_RESET:
            encodeAndSend(PC_OK);
            reset();
            break;
        case PC_MEASURE_ONLY:
            encodeAndSend(PC_OK);
            triggerMeasurements();
            // This contains the state switch to STATE_MEASURE_ADCS.
            break;
        case PC_CHANGE_MEAS_MODE:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_CHANGE_MEASUREMENT_MODE;
            break;
        case PC_STOP:
            currentState = STATE_IDLE;
            encodeAndSend(PC_OK);
            break;
        case PC_SET_VOLTAGE_ONLY:
            takeMeasurements = false;
            waitingForDataFromPC = true;
            initialTime = millis();
            nextVoltageStep = 0;
            currentState = STATE_SET_VOLTAGE;
            break;
        case PC_SET_SERIAL_NR:
            waitingForDataFromPC = true;
            initialTime = millis();
            currentState = STATE_SET_SERIAL_NR;
            break;
    }
    newMessage = false;
}


bool isAllowedCommand(){
    /**
    Check if the received command is among those that the arduino can process.

    Returns
    -------
    True if the message is acceptable
    **/
    // Check that it is one of the understandable commands
    switch(data_received[0]){
        case PC_AUTOGAIN: break;
        case PC_CALIBRATION: break;
        case PC_CONFIGURATION: break;
        case PC_SET_UP_ADCS: break;
        case PC_RESET: break;
        case PC_SET_VOLTAGE: break;
        case PC_MEASURE_ONLY: break;
        case PC_CHANGE_MEAS_MODE: break;
        case PC_STOP: break;
        case PC_SET_VOLTAGE_ONLY: break;
        case PC_SET_SERIAL_NR: break;
        default:
            raise(ERROR_MSG_UNKNOWN);
            return false;
    }
    return true;
}


/** -------------------- PREPARE FOR ENTERING NEW STATES ------------------- **/

void prepareForAutogain(){
    /**Prepare the system to enter a STATE_AUTOGAIN_ADCS.

    Writes
    ------
    adc0Gain, adc1Gain, numMeasurementsDone, summedMeasurements,
    numMeasurementsToDoBackup, numMeasurementsToDo
    **/
    // Reset data, as we will store measurements in the same arrays
    resetMeasurementData();

    // Gains to zero, where we will measure for the gain optimization
    adc0Gain = 0;
    adc1Gain = 0;

    // Set measurement arrays for ripple measurement to zero.
    maximumPeak[0] = 0;
    maximumPeak[1] = 0;
    minimumPeak[0] = 0;
    minimumPeak[1] = 0;

    // Calibrate at the high frequency used for autogain
    selfCalibrateAllADCs(AD7705_500HZ);  // Takes ~13 ms

    // Backup the number of measurement that will be
    // used during normal measurements later on. This
    // is restored at the end of the gain optimization
    numMeasurementsToDoBackup = numMeasurementsToDo;
    numMeasurementsToDo = 25;
}


void triggerMeasurements() {
    /**Trigger the (available) ADCs to start converting.

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel

    Writes
    ------
    adc0Gain, adc1Gain, summedMeasurements, numMeasurementsDone

    Msg to PC
    ---------
    None.

    Goes to state
    -------------
    STATE_MEASURE_ADCS
        At end of execution, effectively starting the measurement
    **/
    // Decrease gain before triggering: need to
    // trigger each time the gain is changed.
    decreaseADCGainsIfNeeded();

    // After triggering, 3/updateRate sec pass
    // before the first data is available
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);

    // Prepare the system for measurement
    resetMeasurementData();
    initialTime = millis();
    currentState = STATE_MEASURE_ADCS;
}


/** ---------------------- STATE & PC-REQUEST HANDLERS --------------------- **/

/** Handler of STATE_GET_CONFIGURATION */
void getConfiguration(){
    /**Send box ID, firmware version, hardware config and serial nr. to PC.
    The serial number is read from the EEPROM.

    Writes
    -----
    hardwareDetected

    Msg to PC
    ---------
    9 data bytes
        the first one is the box ID
        two are the firmware version (M, m)
        two are the hardware configuration as bitmask
        last 4 are the serial number

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_GET_CONFIGURATION
    STATE_IDLE
        Otherwise
    **/
    if (currentState != STATE_GET_CONFIGURATION){
        raise(ERROR_RUNTIME);
        return;
    }

    // Note that the ATmega32U4 of Arduino Micro uses
    // little-endian memory layout, i.e., LSB is
    // at lower memory index

    hardwareDetected.asInt = getHardwarePresent();
    byte serial_nr[4];
    getSerialNumber(serial_nr);
    byte configuration[9] = {BOX_ID,
                             FIRMWARE_VERSION_MAJOR,
                             FIRMWARE_VERSION_MINOR,
                             hardwareDetected.asBytes[1],
                             hardwareDetected.asBytes[0]};
    int address = 0;
    while(address <= 3){
        configuration[address + 5] = serial_nr[address];
        address += 1;
    }
    encodeAndSend(configuration, LENGTH(configuration));
    hardwareNeverChecked = false;
    currentState = STATE_IDLE;
}


// For reasons unknown, the Arduino compiler complains if this
// function is declared after its caller (the next one). Having
// a prototype in the .h file does not help.
void storeAllSelfCalibrationResults(int32_t targetArray[N_MAX_ADCS_ON_PCB][2]) {
    /**
    Store the result of the last self calibration
    of both ADCs (if present) to the array passed

    Parameters
    ----------
    targetArray : 2x2 int32_t
        The array where the calibration results are to be saved.
        The first index identifies which ADC, the second one is
        used to store offset (at index 0) and gain (at index 1).

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel
    **/
    for (int iADC = 0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
        uint16_t adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;      // TODO: using the ADCs struct would make it easier, as you would just need to do externalADCs[iADC].present
        if (hardwareDetected.asInt & adcPresentMask) {
            byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;                    // TODO: similarly, with the ADCs struct, would be enough to do externalADCs[iADC].chipSelect
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: same thing here: externalADCs[iADC].channel
            int32_t offs = AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET);
            int32_t gain = AD7705getCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN);
            targetArray[iADC][0] = offs;
            targetArray[iADC][1] = gain;
        }
    }
}


/** Handler of STATE_CALIBRATE_ADCS */
void calibrateADCsAtAllGains(){
    /**Calibrate (in parallel) the available ADCs for all possible gains.

    The calibration results are stored, so they can be fetched
    later to allow fast switching of the gain.

    Calibration is performed for those update rate and channels
    that come from the PC as the first data message following
    PC_CALIBRATION.

    Each gain value is done over a separate state-loop. This means
    that it takes approximately 360 ms to finish a single run of
    this function (120 ms x 3 points for computing medians), i.e.,
    that the serial line remains unresponsive for ~360 ms.

    Should one need to change the updateRate used during measurements,
    the system must run again through this STATE_CALIBRATE_ADCS.

    Reads
    -----
    msgLength

    Writes
    ------
    initialTime, newMessage, waitingForDataFromPC,
    adc0Channel, adc1Channel, adcUpdateRate,
    calibrationGain, adc0Gain, adc1Gain, selfCalDataVsGain

    Msg to PC
    ---------
    PC_OK.

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_CALIBRATE_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If the ADC channels and update rate data needed from
        the PC do not arrive within 5 s of entering the state
    STATE_ERROR : ERROR_TIMEOUT
        If calibration for a single state-loop iteration (i.e.,
        one value of gain) takes longer than 5 s
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we get too few data bytes, or if either one of the
        channels or the update rate is not an acceptable value
    STATE_CALIBRATE_ADCS (stays)
        Until all the gain values have been calibrated
        for the currently selected channels
    STATE_IDLE
        After calibration is done
    */
    if (currentState != STATE_CALIBRATE_ADCS){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){  // waiting for data from the PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC) {
        // Data has arrived. This part runs only once,
        // the first time the data arrives
        newMessage = false;
        waitingForDataFromPC = false;

        // Check that we got the data that we expected, i.e.,
        // (1) at least: 1 byte for update rate + 1 channel per ADC
        if (msgLength < 1 + N_MAX_ADCS_ON_PCB){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }

        // (2) The update rate is acceptable
        if (not isAcceptableADCUpdateRate(data_received[0])){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }
        adcUpdateRate = data_received[0];

        // (3) The channels are acceptable
        for (byte i = 0; i < N_MAX_ADCS_ON_PCB; i++){
            byte channel = data_received[i+1];
            if (channel != AD7705_CH0 and channel != AD7705_CH1){
                raise(ERROR_MSG_DATA_INVALID);
                return;
            }
        }
        adc0Channel = data_received[1];
        adc1Channel = data_received[2];
    }

    // May timeout in the following if the ADC power is disconnected
    // while we're waiting for a calibration to be finished
    initialTime = millis();

    // The next part runs once every state loop
    if (calibrationGain <= AD7705_MAX_GAIN) {
        // Prepare a place to store values from which medians are
        // calculated: three for median, N_MAX_ADCS_ON_PCB for ADCs,
        // last index is offset(0) & gain(1)
        int32_t selfCalDataForMedian[3][N_MAX_ADCS_ON_PCB][2];

        if (checkIfTimedOut()) {
            adc0Gain = 0;
            adc1Gain = 0;
            return;
        }

        adc0Gain = calibrationGain;
        adc1Gain = calibrationGain;
        // Calibrate 3 times each of the ADCs (in parallel) at
        // the current calibrationGain, and for the adc0Channel
        // and adc1Channel currently selected.
        for (int i = 0; i < 3; i++) {
            selfCalibrateAllADCs(adcUpdateRate);
            storeAllSelfCalibrationResults(selfCalDataForMedian[i]);
            #if DEBUG
                debugMsg("gain=%d, cOffs=%d, cGain=%d",
                         calibrationGain, selfCalDataForMedian[i][0][0],
                         selfCalDataForMedian[i][0][1]);
            #endif
        }

        // Compute the median of the values above, and store
        // it for each ADC, so it can be fetched later
        for (int iADC=0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: using the externalADCs array of struct would make this easier: externalADCs[iADC].channel
            for (int offsetOrGain = 0; offsetOrGain < 2; offsetOrGain++) {
                int32_t median = getMedian(
                    selfCalDataForMedian[0][iADC][offsetOrGain],
                    selfCalDataForMedian[1][iADC][offsetOrGain],
                    selfCalDataForMedian[2][iADC][offsetOrGain]
                    );
                selfCalDataVsGain[calibrationGain][iADC][channel][offsetOrGain] = median;
            }
        }
        calibrationGain++;
    }

    if(calibrationGain > AD7705_MAX_GAIN){
        // Calibration is over.
        // Remember which channels we have calibrated
        for (int iADC=0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
            byte channel = iADC==0 ? adc0Channel : adc1Channel;
            calibratedChannels[iADC][channel] = true;                           // TODO: perhaps should only set to true only for the ADCs that are present. Think what are the consequences.
        }

        // Set gains to the lowest value, and load the
        // new gains and their calibration in the ADCs
        adc0Gain = 0;
        adc1Gain = 0;
        setAllADCgainsAndCalibration();

        // Finally, tell the PC that we're done
        encodeAndSend(PC_OK);
        currentState = STATE_IDLE;
    }
}


/** Handler of STATE_SET_UP_ADCS */
void prepareADCsForMeasurement(){
    /**
    Initialize the ADC from the parameters stored in the first 4 bytes
    of data_received. The bytes have the following meaning:
    - [0 + 1] number of measurements that will be averaged. [0] is the
              high byte, [1] the low byte. Interpreted as an uint16_t.
    - [2] channel of ADC#0
    - [3] channel of ADC#1
    The ADC gains are set to 0 before the calibration is applied.

    Reads
    -----
    data_received

    Writes
    ------
    newMessage, adc0Channel, adc1Channel, numMeasurementsToDo,
    adc0Gain, adc1Gain

    Msg to PC
    ---------
    PC_OK

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_UP_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_UP_ADCS message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive less than the data that we expect (4 bytes)
        or if either channel is not an acceptable value
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_SET_UP_ADCS (stays)
        While waiting for data from the PC
    STATE_IDLE
        Successfully finished
    **/
    if (currentState != STATE_SET_UP_ADCS){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage){  // waiting for data from the PC
        checkIfTimedOut();
        return;
    }

    // Data has arrived
    waitingForDataFromPC = false;
    newMessage = false;

    // Do some checking on the data we received:
    // (1) At least 4 bytes? [2(no. measurements) + 2(channels)]
    if (msgLength < 4){
        raise(ERROR_MSG_DATA_INVALID);
        return;
    }

    numMeasurementsToDo = data_received[0] << 8 | data_received[1];

    // (2) channels acceptable?
    for (byte i = 0; i < N_MAX_ADCS_ON_PCB; i++){
        byte channel = data_received[i+2];
        if (channel != AD7705_CH0 and channel != AD7705_CH1){
            raise(ERROR_MSG_DATA_INVALID);
            return;
        }
    }
    adc0Channel = data_received[2];
    adc1Channel = data_received[3];

    adc0Gain = 0;
    adc1Gain = 0;

    setAllADCgainsAndCalibration();
    if (currentState == STATE_ERROR)   // Some channel was not calibrated
        return;

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}


/** Handler of STATE_SET_VOLTAGE */
void setVoltageWaitAndTrigger(){
    /**
    Ask the DAC to provide a new voltage, if new settings are available.
    Then, wait the prescribed amount of time for the voltage to be stable,
    and trigger the ADCs for a measurement. This may be repeated multiple
    times to allow for multiple voltage steps at once. Triggering occurs
    only after all steps are completed and only if takeMeasurements is
    true.

    For each step, we need 4 bytes from the data_received[] buffer.
    Their meaning is:
    - [0(MSB) + 1(LSB)] DAC value to be set (uint16_t).
    - [2(MSB) + 3(LSB)] dacSettlingTime (uint16_t). Time interval to
            wait before considering the DAC value stable
    Multiple steps are done whenever the data_received[] buffer contains
    multiple 8-bytes consecutive blocks.

    Reads
    -----
    data_received, msgLength

    Writes
    ------
    dacSettlingTime, newMessage, waitingForDataFromPC, summedMeasurements,
    numMeasurementsDone, adc0Gain, adc1Gain, adc0ShouldDecreaseGain,
    adc1ShouldDecreaseGain, initialTime, nextVoltageStep

    Msg to PC
    ---------
    PC_OK, signaling that we're done waiting, and we
    will begin collecting measurements

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_VOLTAGE
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_VOLTAGE message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If we receive data, but it isn't exactly the amount
        we were expecting to get
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If the ADCs should be triggered for measurements but
        one of the ADC channels to be used was not calibrated
    STATE_SET_VOLTAGE (stays)
        While waiting for data from the PC, and until the
        voltage output can be considered stable
    STATE_MEASURE_ADCS
        Successfully finished, accessed via triggerMeasurements()
    STATE_IDLE
        If PC_SET_VOLTAGE_ONLY was received and voltage is set
    **/
    uint16_t dacValue;

    if (currentState != STATE_SET_VOLTAGE){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage and waitingForDataFromPC){
        // Waiting for data from PC
        checkIfTimedOut();
        return;
    }

    if (newMessage and waitingForDataFromPC){ // Do this only once
        // Data has arrived
        waitingForDataFromPC = false;
        newMessage = false;

        // Check that it is the right amount of data
        // For each voltage step: 2 bytes for voltage,
        // 2 bytes for dacSettlingTime
        if (msgLength % 4 != 0){
            raise(ERROR_MSG_DATA_INVALID);

            // And do what we normally would at the
            // beginning of STATE_TRIGGER_ADCS (we're
            // not going there this time)
            resetMeasurementData();
            return;
        }

        dacValue = data_received[0] << 8 | data_received[1];
        dacSettlingTime = data_received[2] << 8 | data_received[3];

        // Set the new voltage
        AD5683setVoltage(CS_DAC, dacValue);
        initialTime = millis();
    }

    // Wait for timer if needed
    if((millis() - initialTime) < dacSettlingTime) return;

    // One voltage step is over.
    nextVoltageStep++;

    // See if we have another one to do
    byte offset = nextVoltageStep*4;
    if (offset < msgLength){
        dacValue = data_received[offset] << 8 | data_received[offset+1];
        dacSettlingTime = data_received[offset+2] << 8 | data_received[offset+3];

        // Set the new voltage
        AD5683setVoltage(CS_DAC, dacValue);
        initialTime = millis();
        return;
    }

    // Finally tell the PC we're done waiting, and
    // trigger the ADCs for measurement if needed
    encodeAndSend(PC_OK);
    if (takeMeasurements){
        triggerMeasurements();  // This switches to STATE_MEASURE_ADCS
    }
    else {
        currentState = STATE_IDLE;
        takeMeasurements = true;
    }
    nextVoltageStep = 0;    // This is not strictly needed, but nicer for cleanup
}


/** Handler of STATE_MEASURE_ADCS */
void measureADCs(){
    /**Acquire measurements from the ADCs.

    (1) Take a defined number of measurements (as per numMeasurementsToDo)
        that are added together in summedMeasurements[].
    (2) When the number of measurement points reaches the requested amount,
        the state of the Arduino is changed to STATE_ADC_VALUES_READY, which
        will handle the communication of the measurement to the PC.

    Reads
    -----
    numMeasurementsDone, numMeasurementsToDo, initialTime

    Writes
    ------
    summedMeasurements, currentState, newMessage

    Msg to PC
    ---------
    None.

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_MEASURE_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If it takes longer than 5s between the beginning of the state
        and completing the number of measurements that need to be averaged
    STATE_ERROR : ERROR_ADC_SATURATED
        If one of the inputs of the ADCs has reached saturation, and
        there is no room to decrease the gain.
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_MEASURE_ADCS (stays)
        Until finished measuring all the data points
    STATE_ADC_VALUES_READY
        Successfully finished
    **/
    if (currentState != STATE_MEASURE_ADCS){
        raise(ERROR_RUNTIME);
        return;
    }

    makeAndSumMeasurements();
    // TODO: if, by chance, the last value measured by one of the ADCs
    //       caused a solid saturation, this is not communicated in any
    //       way to the PC, since the rest of the code goes through, and
    //       the Arduino goes to STATE_ADC_VALUES_READY, swallowing
    //       the error. It could be solved easily by a simple check:
    //       if (currentState == STATE_ERROR) return;
    if (checkIfTimedOut()){
        return;
    }
    if(numMeasurementsDone == numMeasurementsToDo){
        currentState = STATE_ADC_VALUES_READY;
    }
}


/** Handler of STATE_ADC_VALUES_READY */
void sendMeasuredValues(){
    /**Send measurements back to the PC.

    Currently, it actually sends THREE separate messages. Each of
    the messages contains 4 data bytes that, after being re-packed
    (MSB first) correspond to the analog values measured in physical
    units like so:
    - I0(@LEED)   ADC#0, ch0: uAmps
    - HV          ADC#0, ch1: Volts
    - I0(@sample) ADC#1, ch0: uAmps
    - AUX         ADC#1, ch1: Volts
    - LM35 (if present)     : degC
    - LM35 (if not present) : Volts

    The sequence of returned results is: LM35, ADC#0, ADC#1

    Notice that the I0 value (measured at the LEED electronics, i.e.,
    ADC#0, channel 0) is returned in uAmps under the assumption that:
     (i) in 0--2.5V range 1mA produces 1V (i.e., we use 1kOhm)
    (ii) in 0--10V range 1uA produces 1V
    If, a different conversion impedance is used, this correction is
    delegated to the PC.

    Reads
    -----
    fDataOutput
    (indirectly) hardwareDetected, summedMeasurements, millivoltsPerBit,
                 numMeasurementsToDo, adc0Channel, adc1Channel,
                 adc0Gain, adc1Gain

    Writes
    ------
    (indirectly) fDataOutput, summedMeasurements, numMeasurementsDone

    Msg to PC
    ---------
    Three serial messages, each contains 4 bytes (+encoding) with
    values from ADC#0, ADC#1, LM35 (in this order)

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_ADC_VALUES_READY
    STATE_IDLE
        If continuousMeasurement is false
    STATE_MEASURE_ADCS
        If continuousMeasurement is true
    **/
    // TODO: We may later remove this, if we want the option
    // of sending incomplete data back upon request from the PC
    if (currentState != STATE_ADC_VALUES_READY){
        raise(ERROR_RUNTIME);
        return;
    }

    // TODO: Ideally, one would like to rather return a SINGLE MESSAGE
    //       containing 3*4 significant bytes (+ encoding).
    //       In this case, the worst-case scenario message length would be
    //       N_MAX_MEAS*4(bytes)*2(encoding) + 3 (MSG_START, MSG_END, length)
    getFloatMeasurements();

    // Since the ATMega32u4 (Arduino Micro) uses little-endian memory layout, we
    // have flip the bytes over to maintain consistency of our messages, which
    // are in big-endian order
    byte littleToBigEndian[4];
    for (int iADC = 0; iADC < N_MAX_ADCS_ON_PCB+1; iADC++){  // external ADCs + LM35
        for (int i = 0; i < 4; i++){
            littleToBigEndian[i] = fDataOutput[iADC].asBytes[3-i];
        }
        encodeAndSend(littleToBigEndian, 4);
    }
    //encodeAndSend(adc0Gain); //uncomment for debug (adapt python side accordingly)
    resetMeasurementData();
    if (continuousMeasurement) {
        currentState = STATE_MEASURE_ADCS;
        initialTime = millis();
    }
    else {
        currentState = STATE_IDLE;
    }
}


/** Handler of STATE_AUTOGAIN_ADCS */
void findOptimalADCGains(){
    /**Handler of STATE_AUTOGAIN_ADCS.

    Find the optimal gain for all ADCs.

    (1) measure several values. Remain in STATE_AUTOGAIN_ADCS until done.
    (2) for each ADC, pick a gain G such that the worst-case scenario
        measurement is at least ADC_RANGE_THRESHOLD/ADC_POSITIVE_SATURATION
        of the range of values that can be measured with G. The worst-case
        scenario value is calculated from the minimum and maximum values
        measured among the points.

    Reads
    -----
    numMeasurementsToDo, numMeasurementsDone, maximumPeak, minimumPeak,
    hardwareDetected, adc0Gain, adc1Gain

    Writes
    ------
    numMeasurementsDone, adc0Gain, adc1Gain

    Msg to PC
    ---------
    PC_OK : at the state change when receiving the PC_AUTOGAIN command
    adc0Gain, adc1Gain : before returning to STATE_IDLE

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_AUTOGAIN_ADCS
    STATE_ERROR : ERROR_TIMEOUT
        If it takes more than 5s to acquire all the measurements
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_AUTOGAIN_ADCS (stays)
        While it has not yet finished measuring all the values required
    STATE_IDLE
        Successfully finished
    **/
    if (currentState != STATE_AUTOGAIN_ADCS){
        raise(ERROR_RUNTIME);
        return;
    }

    // Notice that we are not triggering the ADCs here, as we are
    // not particularly interested in measuring from a specific
    // point in time. The first data will anyway be available
    // approximately 3/updateRate (i.e., 6 ms @ 500 Hz) after the
    // end of the self-calibration.
    measureADCsRipple();
    if (numMeasurementsDone < numMeasurementsToDo){
        if (checkIfTimedOut()){
            // Do some cleanup:
            resetMeasurementData();
            numMeasurementsToDo = numMeasurementsToDoBackup;

            // Place the ADCs back at gain zero, and
            // restore updateRate and calibration
            adc0Gain = 0;
            adc1Gain = 0;
            setAllADCgainsAndCalibration();
        }
        return;
    }

    // using int32_t because summing int16_t may overflow
    int32_t autogain_value0;
    int32_t autogain_value1;
    autogain_value0 = (max(abs((int32_t)maximumPeak[0]), abs((int32_t)minimumPeak[0]))
                       + abs((int32_t)maximumPeak[0] - (int32_t)minimumPeak[0]));
    adc0RipplePP = maximumPeak[0] - minimumPeak[0];
    autogain_value1 = (max(abs((int32_t)maximumPeak[1]), abs((int32_t)minimumPeak[1]))
                       + abs((int32_t)maximumPeak[1] - (int32_t)minimumPeak[1]));
    adc1RipplePP = maximumPeak[1] - minimumPeak[1];
    // TODO: probably something to check here: if either autogain_value is
    //       already in saturation with gain=0 something is wrong with
    //       the hardware or with the connections. Probably to check only
    //       if the specific ADC is present.
    if(hardwareDetected.asInt & ADC_0_PRESENT){
        while(((autogain_value0 << (adc0Gain + 1)) < ADC_RANGE_THRESHOLD)
                && (adc0Gain < AD7705_MAX_GAIN)){
            adc0Gain++;
        }
    }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        while(((autogain_value1 << (adc1Gain + 1)) < ADC_RANGE_THRESHOLD)
                && (adc1Gain < AD7705_MAX_GAIN)){
            adc1Gain++;
        }
    }

    // Now clean up and update the ADCs with the new gain
    // found, also setting back the normal updateRate, and
    // loading the relevant calibration data.
    resetMeasurementData();
    numMeasurementsToDo = numMeasurementsToDoBackup;
    setAllADCgainsAndCalibration();
    if (currentState == STATE_ERROR)   // Some channel was not calibrated
        return;
    byte ADCgains[2] = {adc0Gain, adc1Gain};
    encodeAndSend(ADCgains, 2);
    currentState = STATE_IDLE;
}


/** Handler of STATE_ERROR */
void handleErrors(){
    /**Clean up after an error, and report it to the PC.

    Should an error occur, the voltage output is always set
    to zero, and the Arduino forgets about the need to
    decrease ADC gains next time.

    Reads
    -----
    currentState, errorTraceback

    Msg to PC
    ---------
    First message : PC_ERROR
    Second message : errorTraceback, i.e., two bytes
        First byte : state of the Arduino while the error occurred
        Second byte : error code

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within a STATE_ERROR
    STATE_IDLE
        Otherwise
    **/
    if (currentState != STATE_ERROR){
        raise(ERROR_RUNTIME);
        return;
    }

    // First, report the error, so the PC knows
    // there may be some cleanup going on
    encodeAndSend(PC_ERROR);
    encodeAndSend(errorTraceback, LENGTH(errorTraceback));

    // Set voltage to zero and forget about gain decrease
    AD5683setVoltage(CS_DAC, 0x0000);
    adc0ShouldDecreaseGain = false;
    adc1ShouldDecreaseGain = false;

    // Then clean up possible mess that caused the error
    switch(errorTraceback[1]){
        case ERROR_SERIAL_OVERFLOW:
            // Discard all characters in the serial buffer,
            // since messages are anyway likely corrupt.
            while(Serial.available()) Serial.read();
            break;
        case ERROR_MSG_TOO_LONG:
            numBytesRead = 0;
            readingFromSerial = false;
            break;
    }
    waitingForDataFromPC = false;
    newMessage = false;
    currentState = STATE_IDLE;
}


/** Handler of PC_RESET */
void reset(){
    /**
    Reset the Arduino, the ADCs and the DACs to the default
    settings, mimicking as much as possible the state right
    after boot-up or a hardware reset. The DAC output is
    also set to zero.

    Writes
    ------
    numBytesRead, msgLength, readingFromSerial, newMessage,
    waitingForDataFromPC, dacSettlingTime, nextVoltageStep,
    hardwareDetected, adcUpdateRate, adc0Channel, adc1Channel
    adc0Gain, adc1Gain, calibrationGain, adc0RipplePP, adc1RipplePP,
    adc0ShouldDecreaseGain, adc1ShouldDecreaseGain,
    numMeasurementsToDo, numMeasurementsToDoBackup,
    numMeasurementsDone, summedMeasurements, calibratedChannels

    Goes to state
    -------------
    STATE_IDLE
    **/
    numBytesRead = 0;
    msgLength = 0;
    readingFromSerial = false;
    newMessage = false;

    currentState = STATE_IDLE;
    waitingForDataFromPC = false;
    continuousMeasurement = false;
    dacSettlingTime = 100;
    nextVoltageStep = 0;

    hardwareDetected.asInt = 0;
    hardwareNeverChecked = true;
    adcUpdateRate = AD7705_50HZ;
    adc0Channel = AD7705_CH0;
    adc1Channel = AD7705_CH0;
    adc0Gain = 0;
    adc1Gain = 0;
    calibrationGain = 0;
    adc0RipplePP = 0;
    adc1RipplePP = 0;
    adc0ShouldDecreaseGain = false;
    adc1ShouldDecreaseGain = false;

    numMeasurementsToDo = 1;
    numMeasurementsToDoBackup = 1;
    resetMeasurementData();
    continuousMeasurementInterval = 0;

    // Rather than clearing the calibration data, we can just
    // reset the flags marking whether channels were calibrated
    // since we always check that they were calibrated before
    // loading calibration data
    for (byte iADC=0; iADC < N_MAX_ADCS_ON_PCB; iADC++){
        calibratedChannels[iADC][0] = false;
        calibratedChannels[iADC][1] = false;
    }

    // Explicitly reset the ADC I/O communication
    AD7705resetCommunication(CS_ADC_0);
    AD7705resetCommunication(CS_ADC_1);

    // Reset the update rate
    AD7705setClock(CS_ADC_0, adcUpdateRate);
    AD7705setClock(CS_ADC_1, adcUpdateRate);

    // As well as channel and gain
    AD7705setGainAndTrigger(CS_ADC_0, adc0Channel, adc0Gain);
    AD7705setGainAndTrigger(CS_ADC_1, adc1Channel, adc1Gain);

    // Reset DAC and set the output voltage to zero
    AD5683reset(CS_DAC);
    AD5683setVoltage(CS_DAC, 0x0000);

    for (int iDevice=0; iDevice<N_MAX_MEAS; iDevice++){
        fDataOutput[iDevice].asFloat = 0.0;
    }
}


/** Handler of STATE_CHANGE_MEASUREMENT_MODE */
void changeMeasurementMode() {
    /**
    Sets the measurement mode either to continous or single measurement.
    To achieve this the continuousMeasurement boolean is set to true or false.
    If the arduino is ordered to do continuous measurements it will set
    numMeasurementsToDo = 1.

    Reads
    -----
    data_received

    Writes
    ------
    continuousMeasurement

    Msg to PC
    ---------
    PC_OK

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_CHANGE_MEASUREMENT_MODE
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_CHANGE_MEAS_MODE message
        and the receipt of data
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If the message is not at least 2 bytes long, or if the received data
        does not describe on or off.
    STATE_CHANGE_MEASUREMENT_MODE (stays)
        While waiting for data from the PC
    STATE_IDLE
        Successfully finished
  **/
    int continuous_mode;
    if (currentState != STATE_CHANGE_MEASUREMENT_MODE){
        raise(ERROR_RUNTIME);
        return;
    }
    if (not newMessage and waitingForDataFromPC){  // waiting for data from the PC
        checkIfTimedOut();
        return;
    }

    // Data has arrived
    waitingForDataFromPC = false;
    newMessage = false;

    if (msgLength < 2){
        raise(ERROR_MSG_DATA_INVALID);
        return;
    }

    continuous_mode = data_received[0];
    // Note that data_received[1] is currently unused
    // but necessary in order to not confuse this data
    // message with a command.

    if (continuous_mode == 1){
        continuousMeasurement = true;
        numMeasurementsToDo = 1;
    }
    if (continuous_mode == 0) {
        continuousMeasurement = false;
    }
    if (continuous_mode != 1 and continuous_mode != 0) {
        raise(ERROR_MSG_DATA_INVALID);
        return;
    }

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}


bool hardwareNotKnown(){
    if(hardwareNeverChecked){
        raise(ERROR_HARDWARE_UNKNOWN);
        return true;
    }
    return false;
}


/** Handler of STATE_SET_SERIAL_NR */
void setSerialNr() {
   /**
    Writes the assigned serial number to the EEPROM.

    Reads
    -----
    data_received

    Msg to PC
    ---------
    PC_OK

    Goes to state
    -------------
    STATE_ERROR : ERROR_RUNTIME
        If this function is not called within STATE_SET_SERIAL_NR
    STATE_ERROR : ERROR_MSG_DATA_INVALID
        If the sent serial number is not 4 bytes long or if the
        sent data contains bytes that do not match the decimal
        ASCII representation of a capital letter or a number.
    STATE_ERROR : ERROR_TIMEOUT
        If more than 5s pass between the PC_SET_SERIAL_NR message
        and the receipt of data.
    STATE_SET_SERIAL_NR (stays)
        While waiting for data from the PC
    STATE_IDLE
        Successfully finished
    **/
    if (currentState != STATE_SET_SERIAL_NR){
        raise(ERROR_RUNTIME);
        return;
    }

    if (not newMessage){  // waiting for data from the PC
        checkIfTimedOut();
        return;
    }

    // Data has arrived
    waitingForDataFromPC = false;
    newMessage = false;

    // Check that we got 4 bytes for the serial number.
    if (msgLength != 4){
        raise(ERROR_MSG_DATA_INVALID);
        return;
    }

    storeSerialNumber(data_received);

    encodeAndSend(PC_OK);
    currentState = STATE_IDLE;
}

/** ---------------- FUNCTIONS CALLED WITHIN STATE HANDLERS ---------------- **/

/** Returns a bit mask of which hardware was detected */
uint16_t getHardwarePresent() {
    int result = 0;
    // Check for ADCs
    // First reset AD7705 IO. Note that (inverted) chip select = high does
    // not reset communication. Reset is done by writing 32 high bits.
    delay(1);     // Make sure that all lines have settled (1 millisec)
    AD7705resetCommunication(CS_ADC_0);
    AD7705resetCommunication(CS_ADC_1);
    delay(1);     // Make sure that all lines have settled (1 millisec)
    byte adc0comm = AD7705readCommRegister(CS_ADC_0, AD7705_CH0);
    if (adc0comm != 0xff) // 0xff if only pullup, no ADC
        result |= ADC_0_PRESENT;
    byte adc1comm = AD7705readCommRegister(CS_ADC_1, AD7705_CH0);
    if (adc1comm != 0xff) // 0xff if only pullup, no ADC
        result |= ADC_1_PRESENT;
    // Check for LM35 temperature sensor: the analog voltage should be within
    // the ADC range and settle to a similar value after connecting to a pullup
    // resistor (the LM35 cannot sink more than about 1 uA, so the pullup will
    // drive it high)
    pinMode(LM35_PIN, INPUT);
    // Unused measurement; the ADC needs time till it works correctly
    analogReadMedian(LM35_PIN);
    delay(10);    // Make sure the voltage has settled (10 millisec)
    int sensorValue0 = analogReadMedian(LM35_PIN);
    // Measure the voltage with internal pullup
    pinMode(LM35_PIN, INPUT_PULLUP);
    delay(10);    // Apply pullup for 10 millisec
    int sensorValue1 = analogReadMedian(LM35_PIN);
    pinMode(LM35_PIN, INPUT);         // Reset for usual measurements
    delay(10);    // Make sure the voltage has settled (10 millisec)
    int sensorValue2 = analogReadMedian(LM35_PIN);
    if (sensorValue0 > 0 && sensorValue0 < LM35_MAX_ADU
            && sensorValue2 > 0 && sensorValue2 < LM35_MAX_ADU
            && abs(sensorValue2 - sensorValue0) < 10
            && sensorValue1 == ARDUINO_ADC_MAX)
        result |= LM35_PRESENT;
    // Check for relay present: If the relay is mounted, it should also have an
    // external pullup that results in about 0.12 V at the pin, about 48 ADUs.
    // If no relay is mounted, it is either open or jumpered to ground,
    // signalling whether the input range is 10 V or 2.5 V range, respectively.
    // If the input is open, the arduino internal pullup would drive it high;
    // if the input is grounded, the voltage should be close to 0.
    // Measure the voltage with arduino pullup
    pinMode(RELAY_PIN, INPUT_PULLUP);
    // Unused measurement just to make sure
    analogReadMedian(RELAY_PIN);
    delay(10);    // Apply pullup for 10 millisec
    int sensorValue3 = analogReadMedian(RELAY_PIN);
    pinMode(RELAY_PIN, INPUT);  // Measure the voltage without pullup
    delay(10);    // No pullup for 10 millisec
    int sensorValue4 = analogReadMedian(RELAY_PIN);
    if (sensorValue3 < ARDUINO_ADC_MAX // Not an open input
            && sensorValue4 > RELAY_MIN_ADU && sensorValue4 < RELAY_MAX_ADU) {
        result |= RELAY_PRESENT;
    }
    else {
        // Check jumper at JP3 indicating 2.5 V I0 range
        // set by user (if no relay)
        pinMode(JP_I0_PIN, INPUT_PULLUP);
        delay(1);
        if (digitalRead(JP_I0_PIN) == 0)
            result |= JP_I0_CLOSED;
        pinMode(JP_I0_PIN, INPUT);     // Pullup off, reduces power consuption
    }
    // Check jumper at JP5 indicating 2.5 V AUX range set by user (if no relay)
    pinMode(JP_AUX_PIN, INPUT_PULLUP);
    delay(1);
    if (digitalRead(JP_AUX_PIN) == 0)
        result |=   JP_AUX_CLOSED;
    pinMode(JP_AUX_PIN, INPUT);        // Pullup off, reduces power consuption
    return result;
}


bool isAcceptableADCUpdateRate(byte updateRate){                                // TODO: may later move to the ADC 'library'
    /**Return whether the given update rate in an acceptable value.*/
    return ((updateRate == AD7705_50HZ)
            or (updateRate == AD7705_60HZ)
            or (updateRate == AD7705_500HZ));
}


void selfCalibrateAllADCs(byte updateRate) {
    /**
    Simultaneous self-calibration for both ADCs
    (if present), using the current channel and gain.

    Parameters
    ----------
    updateRate : {AD7705_50HZ, AD7705_60HZ, AD7705_500HZ}
        The update rate for which calibration has to be performed
        Currently, no check is performed that the update rate
        is actually one of the acceptable values

    Reads
    -----
    adc0Channel, adc0Gain, adc1Channel, adc1Gain
    **/
    if (not isAcceptableADCUpdateRate(updateRate)){
        raise(ERROR_RUNTIME);
        return;
    }

    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705selfCalibrate(CS_ADC_0, adc0Channel, adc0Gain, updateRate);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705selfCalibrate(CS_ADC_1, adc1Channel, adc1Gain, updateRate);
    if (hardwareDetected.asInt & ADC_0_PRESENT)
        AD7705waitForCalibration(CS_ADC_0, adc0Channel);
    if (hardwareDetected.asInt & ADC_1_PRESENT)
        AD7705waitForCalibration(CS_ADC_1, adc1Channel);

    // AD7705 needs about 100--200 us to process the
    // self-calibration result. 1 ms delay is on the
    // safe side, and does not make much of a difference
    // on the overall procedure, which takes ~6/updateRate
    delay(1);
}


void setAllADCgainsAndCalibration() {
    /**
    Set the ADC gains according to the global adc0Gain,
    adc1Gain variables and restores the previously stored
    calibration values from the selfCalDataVsGain array.
    Also makes sure that the update rate of the ADCs is
    the correct one.

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel, adc0Gain,
    adc1Gain, selfCalDataVsGain, adcUpdateRate

    Goes to state
    -------------
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels was not calibrated before
    Same state
        Otherwise
    **/
    for (int iADC = 0; iADC < N_MAX_ADCS_ON_PCB; iADC++) {
        // byte adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;       // BUG: ADC_0_PRESENT/ADC_1_PRESENT are uint16_t, not byte (probably would work with byte, but better have the right type to start with)!
        uint16_t adcPresentMask = iADC==0 ? ADC_0_PRESENT : ADC_1_PRESENT;      // TODO: easier with externalADCs struct: externalADCs[i].present
        if (hardwareDetected.asInt & adcPresentMask) {
            byte chipSelect = iADC==0 ? CS_ADC_0 : CS_ADC_1;                    // TODO: easier with externalADCs[iADC].chipSelect
            byte channel = iADC==0 ? adc0Channel : adc1Channel;                 // TODO: easier with externalADCs[iADC].channel
            byte gain = iADC==0 ? adc0Gain : adc1Gain;                          // TODO: easier with externalADCs[iADC].gain

            // Make sure that the channel has been calibrated before
            if (not calibratedChannels[iADC][channel]){
                raise(ERROR_NEVER_CALIBRATED);
                return;
            }

            int32_t cOffs = selfCalDataVsGain[gain][iADC][channel][0];
            int32_t cGain = selfCalDataVsGain[gain][iADC][channel][1];
            AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_OFFSET, cOffs);
            AD7705setCalibrationRegister(chipSelect, channel, AD7705_REG_GAIN, cGain);

            // Reset the update rate. This is necessary after
            // finishing the auto-gain (that is done at a
            // different frequency), but does not hurt anyway
            AD7705setClock(chipSelect, adcUpdateRate);

            // Trigger to make sure the ADC is up to date with the new
            // calibration data, but we will not read the converted data
            AD7705setGainAndTrigger(chipSelect, channel, gain);
        }
    }
}


// Call this in other measurement functions that we may implement later on!
void decreaseADCGainsIfNeeded(){
    /**
    Decrease the gain of all ADCs if, during the last
    measurement, the value converted by the ADC has
    reached the upper part of the range.

    Writes
    ------
    adc0Gain, adc1Gain, adc0ShouldDecreaseGain, adc1ShouldDecreaseGain

    Goes to state
    -------------
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    Stays the same
        Otherwise
    **/
    if(adc0ShouldDecreaseGain){
        adc0Gain--;
        setAllADCgainsAndCalibration();    // Load the correct calibration
        if (currentState == STATE_ERROR)   // Channel was not calibrated
            return;
        adc0ShouldDecreaseGain = false;
    }
    if(adc1ShouldDecreaseGain){
        adc1Gain--;
        setAllADCgainsAndCalibration();    // Load the correct calibration
        if (currentState == STATE_ERROR)   // Channel was not calibrated
            return;
        adc1ShouldDecreaseGain = false;
    }
}


void makeAndSumMeasurements() {
    /**
    Get one measurement value from all the available ADCs, including the
    LM35, and add them up into summedMeasurements[]. The 'measurements'
    for devices not found are set to 0.

    Whenever a value is measured, it is checked against the saturation
    thresholds, possibly determining a gain switch the next time a
    measurement is performed.

    Notice that this does not trigger the ADCs, that should already be
    ready for measurement before this function is called.
    ADCs have to be triggered after changing the channel.

    Reads
    -----
    hardwareDetected, adc0Channel, adc1Channel

    Writes
    ------
    summedMeasurements, numMeasurementsDone
    (indirectly) adc0Gain, adc1Gain, adc0ShouldDecreaseGain,
                 adc1ShouldDecreaseGain

    Msg to PC
    ---------
    None.

    Goes to state
    -------------
    STATE_ERROR : ERROR_ADC_SATURATED
        If one of the values read by the ADCs reaches a solid
        saturation, and the gain cannot be decreased further
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If one of the ADC channels used was not calibrated
    STATE_MEASURE_ADCS (stays)
        Otherwise
    */
    // TODO: we should discuss what to do when a measurement saturates.
    //   Currently we are adding the measured value anyway. Perhaps we
    //   should shift the summedMeasurements by 1 and not increase the
    //   number of measurements? We could detect this by having
    //   checkMeasurementInADCRange return a boolean (true for in range,
    //   false for out of range), and add the measurement only if it is
    //   in range. In that case, we should take care of skipping completely
    //   this iteration for all the ADCs (i.e. not adding any of the
    //   measurements to summedMeasurements), and not increment the number
    //   of measurements at all (otherwise we screw the averaging).

    // TODO: we are not checking whether this value saturates (to
    //   LM35_MAX_ADU or an appropriate fraction of it, to be discussed),
    //   which may mean that there is a malfunction (e.g., high
    //   temperature). We may want to send an ERROR_TOO_HOT.

    int16_t measurement;
    numMeasurementsDone++;
    if (hardwareDetected.asInt & LM35_PRESENT){
        measurement = analogReadMedian(LM35_PIN);
        summedMeasurements[2] += measurement;
        }   // this one first while we probably have to wait for the others
    if (hardwareDetected.asInt & ADC_0_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
        summedMeasurements[0] += measurement;
        checkMeasurementInADCRange(&adc0Gain, &adc0ShouldDecreaseGain,
                                   measurement, adc0RipplePP);
        }
    if (checkIfTimedOut()){
        return;
    }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        summedMeasurements[1] += measurement;
        checkMeasurementInADCRange(&adc1Gain, &adc1ShouldDecreaseGain,
                                   measurement, adc1RipplePP);
        }
}


void checkMeasurementInADCRange(byte* gain, bool* adcShouldDecreaseGain,
                                int16_t adcValue, int16_t ripple){    // The call to this function would be much easier having externalADCs, as then we would just pass the index of the ADC in the externalADCs array
    /**
    Check whether the (signed) value read is approaching
    the saturation for the current measurement range.

    Parameters
    ----------
    gain : byte
        The current gain of the ADC
    adcShouldDecreaseGain : *bool
        Pointer to the flag of the current ADC, which will be
        used to signal, in the next call to setVoltage(), that
        a gain reduction is necessary.
    adcValue : int16_t
        Signed version of the value measured by the ADC that
        has to be checked for saturation. The check will be
        actually done on abs(adcValue) + abs(ripple).
        (1) If adcValue is in the middle ~50% band of the
            measurement range (i.e., 2*ADC_RANGE_THRESHOLD),
            the value is still acceptable, but gain will soon
            need to be decreased. Set adcShouldDecreaseGain.
        (2) If adcValue is in solid saturation (full scale),                    // TODO: can't we avoid redoing all measurements? Perhaps it's enough to shift the correct summedMeasurements[i] by >> 1?
            the measurement itself is incorrect, and needs to
            be repeated from scratch after decreasing the gain.
        (3) If adcValue is in solid saturation but the gain
            cannot be decreased, go to STATE_ERROR,
            with ERROR_ADC_SATURATED code.
    ripple : int16_t
        Peak-to-peak ripple voltage measured at the lowest
        gain for the ADC in question. This value is the one
        set while in STATE_AUTOGAIN_ADCS.

    Goes to state
    -------------
    STATE_ERROR : ERROR_ADC_SATURATED
        If adcValue is in solid saturation, and the gain cannot
        be decreased. Otherwise, the currentState is preserved.
    STATE_ERROR : ERROR_NEVER_CALIBRATED
        If the ADC channel used was not calibrated
    Stays unchanged
        Otherwise
    **/
    if (abs(adcValue) > (ADC_RANGE_THRESHOLD - abs(ripple>>*gain))
            && (*gain > 0)
            && !(*adcShouldDecreaseGain)){
        // The measured value is above the "saturation" threshold,
        // but not yet at true saturation, which would make the
        // measured value completely wrong. Defer the decrease of
        // gain to someone else: currently this is done only when
        // setVoltage() gets called, i.e. at the next energy step.

        // NB: using the absolute value of the signed measurement
        //     means that we deem the value 'above threshold'
        //     when its unsigned version is outside a band twice
        //     as large as the one used in findOptimalADCGains()
        *adcShouldDecreaseGain = true;
    }

    if(abs(adcValue) >= ADC_SATURATION){
        if(*gain > 0){
            (*gain)--;
            setAllADCgainsAndCalibration();
            if (currentState == STATE_ERROR)   // Some channel was not calibrated
                return;
            resetMeasurementData();
            *adcShouldDecreaseGain = false;
        }
        else{
            raise(ERROR_ADC_SATURATED);
        }
    }
}


void resetMeasurementData() {
    /** Reset summedMeasurements and numMeasurementsDone to 0 **/
    for (int i = 0; i < LENGTH(summedMeasurements); i++)
        summedMeasurements[i] = 0;
    numMeasurementsDone = 0;
}


void getFloatMeasurements() {
    /**
    Convert the global summedMeasurements to float in physical units
    (Volts, uAmps, degC), setting the value in the fDataOutput array.

    Reads
    -----
    hardwareDetected, summedMeasurements, millivoltsPerBit, numMeasurementsToDo,
    adc0Channel, adc0Gain, adc1Channel, adc1Gain

    Writes
    ------
    fDataOutput
    **/
    if (hardwareDetected.asInt & ADC_0_PRESENT) {
        fDataOutput[0].asFloat = (
            summedMeasurements[0] * millivoltsPerBit[adc0Gain]
            / numMeasurementsDone
            );
        // ADC#0, channel 0: I0 input (returned as milliVolts).
        // The actual voltage at the input can be in either 0-2.5 V
        // (jumper closed) or 0-10 V (jumper open) ranges. This means
        // scaling the value by ADC_0_CH0_SCALE_JO if the jumper is open,
        // and doing nothing if it's closed (ADC range is already 0-2.5V)
        if (adc0Channel == AD7705_CH0) {
            if ((hardwareDetected.asInt & JP_I0_CLOSED) == 0)  // 0-10 V with jumper open
                fDataOutput[0].asFloat *= ADC_0_CH0_SCALE_JO;
        }
        // ADC#0, channel 1: high voltage (Volts)
        else
            fDataOutput[0].asFloat *= ADC_0_CH1_SCALE;
    }

    if (hardwareDetected.asInt & ADC_1_PRESENT) {
        fDataOutput[1].asFloat = (
            summedMeasurements[1] * millivoltsPerBit[adc1Gain]
            / numMeasurementsDone
            );
        // ADC#1, channel 0: I0 at biased sample (microAmps)
        if (adc1Channel == AD7705_CH0) {
            fDataOutput[1].asFloat *= ADC_1_CH0_SCALE;
        }
        else{
            // ADC#1 channel 1: AUX
            // The actual voltage at the input can be in either 0-2.5 V
            // (jumper closed) or 0-10 V (jumper open) ranges. This means
            // scaling the vale by JP_AUX_CLOSED if the jumper is open,
            // and doing nothing if it's closed (ADC range is already 0-2.5V)
            if ((hardwareDetected.asInt & JP_AUX_CLOSED) == 0)  // 0-10 V with jumper open
                fDataOutput[1].asFloat *= ADC_1_CH1_SCALE_JO;
        }
    }

    if (hardwareDetected.asInt & LM35_PRESENT) // LM35: degrees C
        fDataOutput[2].asFloat = summedMeasurements[2] * ARDUINO_ADU_DEG_C / numMeasurementsDone;
}


void measureADCsRipple(){
    /**
    Measure a single value from any of the ADCs present, and
    update the maximum and minimum values measured so far.

    Reads
    -----
    hardwareDetected

    Writes
    ------
    maximumPeak, minimumPeak, numMeasurementsDone
    **/
    int16_t measurement;
    if (hardwareDetected.asInt & ADC_0_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_0, adc0Channel);
        if(measurement > maximumPeak[0]){
            maximumPeak[0] = measurement;
            }
        if(measurement < minimumPeak[0]){
            minimumPeak[0] = measurement;
            }
        }
    if (hardwareDetected.asInt & ADC_1_PRESENT){
        measurement = AD7705waitAndReadData(CS_ADC_1, adc1Channel);
        if(measurement > maximumPeak[1]){
            maximumPeak[1] = measurement;
            }
        if(measurement < minimumPeak[1]){
            minimumPeak[1] = measurement;
            }
        }
    numMeasurementsDone++;
}
