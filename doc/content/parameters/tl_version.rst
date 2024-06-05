.. _tl_version:

TL_VERSION
==========

TL_VERSION specifies which version of the TensErLEED code should be used. 
ViPErLEED will look for TensErLEED code in the ``tensorleed/TensErLEED-vX.XX`` 
directory, where X.XX is the version number (interpreted as float).

**Default:** Picks the highest version number found in the ``tensorleed`` 
directory.

**Allowed values:** positive float

**Syntax:**

::

   TL_VERSION = 1.6

If the specified version is not found, execution will stop.
