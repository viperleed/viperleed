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

Version numbers use the canonical format ``v<major>.<minor>.<patch>`` to
denote major, minor and patch updates. The ``v`` prefix is optional.
Versions prior to ``v2.0.0`` will also be recognized when given in the
legacy format of ``v<major>.<minor><patch>``.

.. versionchanged:: 0.14.2
    Added support for TensErLEED v2.1.0. TensErLEED versions v1.7.2 to v1.7.6
    and v2.0.0 compute wrong tensors for the topmost layer, which affects
    delta amplitudes, structure searches and error curves (but not the I(V)
    curves of reference calculations). A warning is now logged when one of
    these versions is used. We recommend upgrading to TensErLEED v2.1.0.
