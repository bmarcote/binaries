import numpy as np
from astropy import units as u
from astropy import time
from astropy import coordinates as coord




def sky_offset_from_proper_motion(source_coord, epoch, ref_epoch=None):
    """Given a source position (which can have a defined proper motion inside),
    it determines the position at the given epoch due to the proper motion.

    Parameters:
        source_coord : astropy.coordinates.SkyCoord
            The sky coordinates of the source to evaluate. It must include proper motion info.
        epoch : astropy.time.Time
            Epoch at which the offset position must be determined.
        ref_epoch : astropy.time.Time [OPTIONAL]
            Epoch at which the offsets are referred to. If None, the obstime in source_coord
            will be taken as reference.
    Returns
        offset : (ra_offset: float, dec_offset: float)
            Tuple with the expected offsets in RA and DEC (astropy.units are provided).
    """
    delta_ra = source_coord.ra + source_coord.pm_ra_cosdec*(epoch - source_coord.obstime)
    delta_dec  = source_coord.dec + source_coord.pm_dec*(epoch - source_coord.obstime)
    if ref_epoch is not None:
        delta_ra0 = source_coord.ra + source_coord.pm_ra_cosdec*(ref_epoch - source_coord.obstime)
        delta_dec0  = source_coord.dec + source_coord.pm_dec*(ref_epoch - source_coord.obstime)
        return delta_ra - delta_ra0, delta_dec - delta_dec0

    return delta_ra, delta_dec



def sky_offset_from_parallax(source_coord, epoch, parallax=None):
    """Given a source position, determines the position offset in the sky
    due to the parallax motion at the given epoch.

    This code is based on a Fortran-77 script (fit_parallax_multi_4d.f) from
    Katharina Immer.

    Parameters:
        source_coord : astropy.coordinates.SkyCoord
            The sky coordinates of the source to evaluate.
        epoch : astropy.time.Time
            Epoch at which the offset position must be determined.
        parallax : unit.Quantity or float [OPTIONAL]
            Parallax of the source. Required only if distance is not provided
            in source_coord or needs to be ignored. If units are not provided,
            milliarcsecond are assumed.

    Returns
        offset : (ra_offset: float, dec_offset: float)
            Tuple with the expected offsets in RA and DEC (astropy.units are provided).
    """
    obliquity = 23.439*u.deg
    # Number of dates from 2000.0
    n_days = epoch.jd - 2451545.0
    # Mean longitude of the Sun (corrected for aberration of light)
    long_sun = (280.460 + 0.9856474*n_days)*u.deg
    # Cartesian coordinates of the Sun at the given epoch
    x = np.cos(long_sun)
    y = np.sin(long_sun)*np.cos(obliquity)
    z = np.sin(long_sun)*np.sin(obliquity)
    # Correcting for eccentricity of the Earth orbit
    long_earth =  2*np.pi*(epoch.byear - 0.257)
    factor = 1.0 + 0.0167*np.sin(long_earth)
    proj_ra = y*np.cos(source_coord.ra) - x*np.sin(source_coord.ra)
    proj_dec = z*np.cos(source_coord.dec) - x*np.cos(source_coord.ra)*np.sin(source_coord.dec)\
               - y*np.sin(source_coord.ra)*np.cos(source_coord.dec)
    if parallax is not None:
        if not isinstance(parallax, u.Quantity):
            parallax = parallax*u.mas
    else:
        if source_coord.distance.unit is u.dimensionless_unscaled:
            raise u.UnitConversionError('Distance is not set and thus parallax cannot be determined.')
        parallax = source_coord.distance.to(u.mas, equivalencies=u.parallax())

    return (factor*parallax*proj_ra, factor*parallax*proj_dec)



# Functions to determine orbital apparent motions
def mean_anomaly(source, epoch):
    """Returns the mean anomaly ot a given epoch.
    Parameters:
        source : Sources.binaries.System
            The binary source to compute the motion. The following orbital parameters must
            be provided:
                - period : orbital period of the system.
                - t0 : reference epoch at which origin of phases take place.
                - periastron_phase : the orbital phase at which periastron takes place.
        epoch : astropy.time.Time
            Epoch at which the mean anomaly must be determined.
    Returns:
        mean_anomaly : float or astropy.units.Quantity
            The mean anomaly at the given epoch.
    """
    return np.fmod(((epoch.mjd-source.t0.value)*u.day/source.period.value).decompose()-source.periastron_phase.value, 1.0)


def eccentric_anomaly(source, mean_anomaly):
    """Returns the eccentric anomaly ot a given epoch. Computes it by a simple iteraction process.
    Parameters:
        source : Sources.binaries.System
            The binary source to compute the motion. The following orbital parameters must
            be provided:
                - eccentricity : eccentricity of the orbit
        mean_anomaly : float or astropy.units.Quantity
            The mean anomaly to use.
    Returns:
        eccentric_anomaly : float or astropy.units.Quantity
            The eccentric anomaly at the given mean anomaly.
    """
    ecc_anom0 = mean_anomaly*2*np.pi
    ecc_anom = 2*np.pi*mean_anomaly + source.eccentricity.value*np.sin(ecc_anom0*u.rad)
    while (np.maximum.reduce(np.fabs(ecc_anom - ecc_anom0)) > 5e-8):
        ecc_anom0 = ecc_anom.copy()
    ecc_anom = 2*np.pi*mean_anomaly + source.eccentricity.value*np.sin(ecc_anom0*u.rad)
    return np.remainder(ecc_anom/(2*np.pi), 1.0)


def true_anomaly(source, eccentric_anomaly):
    """Returns the true anomaly ot a given eccentric anomaly.

    Parameters:
        source : Sources.binaries.System
            The binary source to compute the motion. The following orbital parameters must
            be provided:
                - eccentricity : eccentricity of the orbit
        eccentric_anomaly : float or astropy.units.Quantity
            The eccentric anomaly at the given mean anomaly.
    Returns:
        true_anomaly : float or astropy.units.Quantity
            The true anomaly at the given eccentric anomaly.
    """
    true_anom = (np.arctan(np.sqrt((1+source.eccentricity.value)/(1-source.eccentricity.value))*np.tan(eccentric_anomaly*2*np.pi*u.rad/2.))/np.pi).value
    if np.isscalar(true_anom):
        return true_anom + 1.0 if true_anom < 0 else true_anom

    true_anom[np.where(true_anom < 0.0)] += 1.0
    return true_anom


def apparent_separation(source, eccentric_anomaly):
    """Returns the apparent separation in the sky between the focus of the orbit and the
    position of the object at the given eccentric anomaly.

    Parameters:
        source : Sources.binaries.System
            The binary source to compute the motion. The following orbital parameters must
            be provided:
                - eccentricity : eccentricity of the orbit
                - inclination : inclination of the orbit with respect to the plane of the sky.
                - semimajor_axis : semimajor axis.
                - omega : argument of periastron (angle from the ascending node to periastron).
        eccentric_anomaly : float or astropy.units.Quantity
            The eccentric anomaly at which the apparent separation will be computed.
    Returns:
        apparent_separation : float or astropy.units.Quantity
            The apparent separation of the source at the given eccentric anomaly.
    """
    radius = source.semimajor_axis.value*(1 - source.eccentricity.value*np.cos(eccentric_anomaly*2*np.pi*u.rad))
    true_anom = true_anomaly(source, eccentric_anomaly)*2*np.pi*u.rad
    return radius*np.sqrt((np.sin(true_anom + source.omega.value)*np.cos(source.inclination.value))**2 + \
                          np.cos(true_anom + source.omega.value)**2)


def apparent_position_angle(source, true_anomaly):
    """Returns the apparent position angle in the sky between the focus of the orbit and the
    position of the object at the given true anomaly.

    Parameters:
        source : Sources.binaries.System
            The binary source to compute the motion. The following orbital parameters must
            be provided:
                - inclination : inclination of the orbit with respect to the plane of the sky.
                - omega : argument of periastron (angle from the ascending node to periastron).
                - Omega : longitude of the ascending node.
        true_anomaly : float or astropy.units.Quantity
            The true anomaly at which the apparent separation will be computed.
    Returns:
        apparent_position_angle : float or astropy.units.Quantity
            The apparent position angle of the source at the given true anomaly.
    """
    pa = np.arctan2(np.sin(true_anomaly*2*np.pi*u.rad+source.omega.value)*np.cos(source.inclination.value),\
            np.cos(true_anomaly*2*np.pi*u.rad+source.omega.value)) + source.Omega.value
    if np.isscalar(pa.value):
        return pa + 2*np.pi*u.rad if pa < 0.0 else pa

    pa[np.where(pa < 0.0)] += 2*np.pi*u.rad
    return pa


def apparent_delta_ra_dec(rho, theta):
    """Converts an offset given as (rho, theta) to the one in (delta_RA, delta_DEC).

    Parameters:
        rho : float or astropy.units.Quantity
            The apparent radial separation of the source.
        theta : float or astropy.units.Quantity
            The apparent position angle of the source.
    Returns:
        offset : (ra_offset, dec_offset)
            The position offset of the object in RA and DEC.
    """
    return (rho*np.sin(theta), rho*np.cos(theta))


def sky_offset_from_orbit(source, epoch):
    """Given a binary source, determines the position offset in the sky at the
    given epoch due to the orbital motion.

    Parameters:
        source : Sources.binaries.System
            The binary source to compute the motion. The following orbital parameters must
            be provided:
                - eccentricity : eccentricity of the orbit
                - inclination : inclination of the orbit with respect to the plane of the sky.
                - semimajor_axis : semimajor axis.
                - t0 : reference epoch at which origin of phases take place.
                - periastron_phase : the orbital phase at which periastron takes place.
                - omega : argument of periastron (angle from the ascending node to periastron).
                - Omega : longitude of the ascending node.
                - period : orbital period of the system.

        epoch : astropy.time.Time
            Epoch at which the offset position must be determined.

    Returns:
        offset : (ra_offset: float, dec_offset: float)
            Tuple with the expected offsets in RA and DEC (astropy.units are provided).
    """
    # Equations obtained from Descomps (2005).
    ecc_anom = eccentric_anomaly(source, mean_anomaly(source, epoch))
    true_anom = true_anomaly(source, ecc_anom)
    # To angular scales
    rho = (apparent_separation(source, ecc_anom).to(u.AU).value/source.distance.value.to(u.kpc).value)*u.mas
    theta = apparent_position_angle(source, true_anom)
    # Get the RA and DEC coordinates
    return apparent_delta_ra_dec(rho, theta)




