
.. _developer:

Reference Files
---------------
The HST reference files in config/hst may occasionally need to be remade with
updated values from the instrument teams. 
A convenience generator for the reference files is provided in 
`astrogrism.HST.create_HST_reffiles.py`, which creates the specws and wavelengthrange
reference files based on a GRISMCONF configuration file. The distortion generator 
is currently in progress. To generate the files, provide the path to the GRISMCONF 
configuration file and the name of the grism as arguments to the 
`create_reference_files` function:

    from astrogrism.config.HST import create_reference_files
    create_reference_files("G102.conf", 'G102')