.. _tl_version:

TL_VERSION
==========

TL_VERSION specifies which version of the TensErLEED code should be used.
ViPErLEED will look for TensErLEED code in the directory specified by the
``--tensorleed`` argument or the ``VIPERLEED_TENSORLEED`` environment variable
(see the :ref:`How to run section<how_to_run>`). If ``TL_VERSION`` is not set,
ViPErLEED will use the latest version found in the specified directory.

**Default:** Pick the latest version number found in the ``tensorleed``
directory.

**Allowed values:** Valid version numbers (e.g., ``2.0.0``)

**Syntax:**

::

   TL_VERSION = 2.0.0
   TL_VERSION = v1.7.6

If the specified version is not found, execution will stop.
