.. include:: /substitutions.rst

.. _tl_ignore_checksum:

TL_IGNORE_CHECKSUM
==================

TL_IGNORE_CHECKSUM toggles SHA-256 checksum calculation for TensErLEED source
code files. During |calc| execution, various parts of TensErLEED are compiled
at run-time. If these source code files were altered, this may lead to
unintended behavior at best and security vulnerabilities at worst. To safeguard
against this (at least in part), |calc| can automatically perform checksum
calculations on the source code files and compare them with known checksums
before compilation. If TL_IGNORE_CHECKSUM is set to True, this will be skipped.
Performing the automated checksums is recommended, but be aware that this not a
completely fail-safe solution. If you are relying on it, ensure that the
checksums stored in |calc| (calc/lib/_checksums.dat) were not tampered with.

Users may want to disable checksum validation for the purposes of custom
modifications to TensErLEED code or for high-throughput calculations.

**Default:** False

**Allowed values:** True/False

.. admonition:: Syntax

   ::

      TL_IGNORE_CHECKSUM = True

Changelog
---------

.. versionchanged:: 0.11.0
   Changed the default value from True to False.
