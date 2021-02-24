# viperleed

Documentation of the API and how to use the package can be found at: https://www.iap.tuwien.ac.at/www/protected/surface/leediv/index

If you do not have access, please contact michele.riva@tuwien.ac.at

This repository includes the following main executables:
- GUI: offers a LEED pattern preview, and exports a "pattern file" required by the spot tracker
- tleedm: the "TensErLEED Manager". Performs LEED-IV calculations based on the TensErLEED package from input files as documented in the wiki. This creates quite a lot of files, so it is recommended to have an additional job script (see e.g. example-job.sh) moving the input files and tleedm executable into a work folder and executing it there. When exiting, tleedm will create a "manifest" file, which is a list of files and folders that should be copied back.
- bookkeeper: a helper utility that copies output from a previous tleedm run into a directory "history", and collects information about previous runs in a file "history.info". See example-job.sh for an example of bookkeeper usage. Execute bookkeeper with flag "-c" for continuation jobs, i.e. to overwrite the POSCAR and VIBROCC input files in the main folder with the latest POSCAR_OUT and VIBROCC_OUT from the OUT folder.
