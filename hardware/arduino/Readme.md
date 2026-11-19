# ViPErLEED Hardware -- Control firmware

Here you find the source code that is used as a 'firmware' for the Arduino controller. The source is contained in folder './viper-ino', which is structured so that the code can be uploaded via the Arduino IDE.
The code is split into the main file (viper-ino.h, viper-ino.h) and into two 'libraries' for the DAC and ADC chips that we are currently using (as of hardware v8).

The logics used in the Arduino firmware, including all commands, as well as errors (and their codes) can be found in './doc/State_Machine.pdf' (the './doc' subfolder also contains the .tex source used to generate the documentation). **Should anyone ever change the logics and/or add new functionality it is very advisable to update the documentation right away.**

As an alternative to compilation and upload via the IDE, the python module 'upload-sketch.py' provides basic functionality to
 * Download the Arduino Command-Line Interface (CLI) from GitHub, picking the version more appropriate for your operating system
 * Extract the archive to a folder './arduino-cli'
 * Compile, check and upload the code to the Arduino Micro board
 
At the present stage (2021-06-10), the python module should be considered as a stub, as it needs significant polishing and testing.

Currently (2021-06-10), the folder './temp-python-control' contains the terminal-based interface for controlling the ViPErLEED hardware and a Camera. It will eventually be merged within the viperleed.guilib module.
