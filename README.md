# binaries

This project is intended to simplify the calculus of parameters and ephemerides from binary systems, like the epoch at which a particular orbital phase takes place or the determination of the orbit in the plane of the sky given the orbital parameters of a system.


## Orbital and system parameters

All the information to be used in the program and related to each particular binary system will be stored as a `*.ini` file inside the `systems` folder. See one of the current examples to create files for new systems.

For each system (aka binary) you can define as many parameters as you wish. all of them will be recognized when loading for that particular system. Some of them must be required depending on different functions to be executed. Note that RA and DEC are the only mandatory parameters to initialize binaries.


## Main program: binaries.py


Reads all the binaries available under the `systems` directory and stores them inside a `Binaries` object. You can access to a particular system with `Binaries().get_binary(binary_name)`, which will contain the coordinates as an `astropy.SkyCoord` parameter, and all the other parameters defined in the ini file.



# To Do

- Modify `binaries.py` to support proper motions from the init files (the `SkyCoord` part must be modified.
  If it is provided, then a reference epoch should be provided.
- Implement `astropy.time` in `date2phase.py` and `binaries.py`. It natively supports MJD conversions.


