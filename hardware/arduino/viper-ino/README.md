### Development notes
1. Changing communication markers must be accompanied by an increase of
	the **major** of the firmware version.
2. Adding commands that change the current behavior of the
	ViPErinoController must be accompanied by an increase of the
	**major** of the firmware version. (e.g., hardware detection
	or changes to the state machine)
3. New commands that do not change the behavior and do not affect
	communication protocol do not require a change of the firmware
	**major**. They must be must be accompanied by an increase of
	the **minor** of the firmware version. (e.g., supporting the relay)
4. Minor bug fixes must be accompanied by an increase of the **minor**
	of the firmware version.
