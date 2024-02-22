.. _tl_ignore_checksum:

TL_IGNORE_CHECKSUM
==================

TL_IGNORE_CHECKSUM toggles SHA-256 checksum calculation for TensErLEED source code files. During tleedm execution, various parts of TensErLEED are compiled at run-time. If these source code files were altered, this may lead to unintended behavior at best and security vulnerabilities at worst. To safeguard against this (at least in part), tleedm can automatically perform checksum calculations on the source code files and compare them with known checksums before compilation. If TL_IGNORE_CHECKSUM is set to True, this will be skipped. Performing the automated checksums is recommended, but be aware that this not a completely fail-safe solution. If you are relying on it, ensure that the checksums stored in tleedm (tleedmlib/files/checksums.py) were not tampered with.

Users may want to disable checksum validation for the purposes of custom modifications to TensErLEED code or for high-throughput calculations.

**Default:** True (for testing, will be set to False in a future release)

**Allowed values:** True/False

**Syntax:**

::

   TL_IGNORE_CHECKSUM = False
