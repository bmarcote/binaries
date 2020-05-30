# Binaries

This project is intended to simplify the calculus of parameters and ephemerides from binary systems, like the epoch at which a particular orbital phase takes place or the determination of the orbit in the plane of the sky given the orbital parameters of a system.

Currently the project is just a collection of independent scripts focused on different utilities. However, the code is under heavy development and the idea is to obtain a consistent and unified package.  

This project is currently composed of different scripts that can work in an isolated way:



## MJD

This module (`mjd.py`) allows the user a quick conversion between different time units widely used in astrophysics: Julian Day (JD), Modified Julian Day (MJD), and gregorian time (Python `datetime` objects).

The program has different functions to convert between one and another, like `date2mjd` or `mjd2jd`.



## Binaries

The module `binaries.py` retrieves the information stored for different binaries under the `data` directory (currently contains information from well-known gamma-ray binaries and a few X-ray binaries) in a `config` (`*.ini`) file.

The user add as many binaries as wished, and can define as many parameters as wished, all of them will be recognized when loading for that particular system. Note that only some of them may be required depending on the different functions to be executed. The coordinates of the source  (`ra`, `dec`) are the only mandatory parameters to initialize binaries.

One can load the `binaries.Binary()` to create an object containing all binaries in the `data` directory. Then one can get the properties of one specific binary by calling `.get_binary({*binary_name*})`.



## Date to phase

The `date2phase.py` module allows the user a quick conversion between dates and orbital phases for defined binaries (`binaries.Binary` objects). 



## Colliding Wind Binaries

The module `colliding_wind_binaries.py` allows the user to determine the contact discontinuity (CD) curve expected from the shock front produced by the wind collision of the two stars in the binary. We followed the analytical form from [Canto et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...469..729C/abstract).

This module allows the user to extract the relevant physical information from this binaries. One can user the traditional description of the CD with the class `WCR`, or compute the sky coordinates where the CD should be placed with `WCR_rs`. 



## Read before using the code

Note that this project combines different independent scripts. An unification of the code is yet to be done. The `binaries` part is an old code that can be highly improved and modernized in the coming future. The `date2phase.py` code can take advantage of newer packages as `astropy` (not developed at the time).



## Dependencies

This program has been written in `Python 3` and requires the following dependencies:

- `numpy`.
- `scipy`.
- `astropy`.



