.. _log_search:

LOG_SEARCH
==========

LOG_SEARCH defines whether the standard output from the TensErLEED search 
program should be written as search-``[timestamp]``.log, or whether the 
output should be routed to /dev/null. The search log normally contains 
nothing of value, so this is mostly useful for debugging purposes.

**Default**: LOG_SEARCH = True (write log)

**Syntax**:

::

   LOG_SEARCH = False

**Acceptable values**: True, False, T, F (not case sensitive)

Note that for TensErLEED versions 1.6 and below, the search log file can be 
very large. In this case, even if LOG_SEARCH = True, the log file will not 
be added to the manifest.
