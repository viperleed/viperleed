__authors__ = (
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-07-15"
__license__ = "GPLv3+"

import multiprocessing as mp

# Set the multiprocessing start method to 'spawn' to ensure compatibility
# with the JAX backend and to avoid issues with forked processes.
# This should not affect other parts of the code, but may cause issues when not
# set and using JAX.

def setup_multiprocessing_method_spawn(force=False):
    if not mp.parent_process():  # Only in the main process
        mp.set_start_method("spawn", force=force)