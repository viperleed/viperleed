.. job_script:

The job script
==============

The job script, typically either ``job.py`` or ``job.sh`` is the entry point for a ViPErLEED run.
Examplary job scripts are provided with ViPErLEED that need minimal adjustments.

To run a normal ViPErLEED calculation, copy the provided job script into the calculation's source directory and fill in the variables ``vpr_path`` and ``work_path`` with the filepaths to the directory containing the ViPErLEED source code and the desired location of the ``work`` directory, respectively. See also :ref:`How To Run<how_to_run>`.

With the :term:`Python` version of the job script (``job.py``), you can also specify the source and work directories with the ``-s`` and ``-w`` flags respectively.

**TODO (AMI):** I think it would be good to reorganize tleedm such that the job_script becomes more or less redundant: move most stuff from job.py to the main() of tleedm.py and make tleedm a callable command. Also, while we are at it, maybe incorporate the bookkeeper into tleedm/job too. I don't see why one wouldn't want to run it (possibly make it a parameter).