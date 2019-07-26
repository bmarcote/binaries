"""Convertions between dates and orbital phase for binary systems.
"""

import numpy as np
import datetime as dt
from astropy import units as u
import mjd
import binaries

def get_system_parameters(binary_system, superorbital=False):
    """Returns the orbital parameters (T0, P and p.periastro) of the given
    binary systems.

    Arguments
        binary_system : binaries.System
            The binary system to be considered.
        superorbital : bool [default: False]
            In case the system has a superorbital period, considers it instead
            of the orbital one.

    Returns
        T0 : float
            The epoch of reference for the orbit.
        Period : float
            The orbital period of the orbit.
        phi0 : float
            The orbital phase at which the periastron takes place.

    Raises
        KeyError -- in case any of the parameters is not found in the system
            (it searches for t0, period, periastron_phase.
    """
    # Copy the parameters from the imported file
    t0, p, phi = None, None, None
    try:
        if superorbital:
            t0 = binary_system.t0_superorbital.value
            p = binary_system.period_superorbital.value.to(u.d).value
            phi = binary_system.periastron_phase.value
        else:
            t0 = binary_system.t0.value
            p = binary_system.period.value.to(u.d).value
            phi = 0.
        return t0, p, phi

    except KeyError:
        raise KeyError('ERROR: Attributes t0, period and/or periastron_phase not found.')


def date2phase(epoch, binary_system, superorbital=False):
    """Return the orbital and post-periastron phase of the given binary system
    for a given epoch entered as datetime.

    Arguments
        epoch : datetime
            Epoch at which the orbital phases will be computed.
        binary_system : binaries.System
            The binary system to be considered.
        superorbital : bool [default: False]
            In case the superorbital period (instead of the normal one) need to be computed.

    Returns
        phase : float
            The orbital phase at which the system is at the given epoch.
        phase_postperiastron : float
            The phase post-periastron at which the system is at the given epoch.
    """
    return mjd2phase(mjd.date2mjd(epoch), binary_system, superorbital=superorbital)


def mjd2phase(epoch, binary_system, superorbital=False):
    """Return the orbital and post-periastron phase of the given binary system
    for a given epoch entered as Modified Julian Date.

    Arguments
        epoch : datetime
            Epoch at which the orbital phases will be computed.
        binary_system : binaries.System
            The binary system to be considered.
        superorbital : bool [default: False]
            In case the superorbital period (instead of the normal one) need to be computed.

    Returns
        phase : float
            The orbital phase at which the system is at the given epoch.
        phase_postperiastron : float
            The phase post-periastron at which the system is at the given epoch.
    """
    t0, p, per = get_system_parameters(binary_system, superorbital=superorbital)
    phase = (epoch - t0)/p % 1
    phase_postperiastron = phase
    if per != 0 and per != None:
        phase_postperiastron = (phase - per + 1) % 1
    return phase, phase_postperiastron


def phase2dates(phase, binary_system, initdate, enddate, postperiastron=False,
                superorbital=False):
    """Find all the dates at which the binary system is in a determinate phase
    for a given range of dates. You can choose if that phase is or not a phase
    refered at the periastron (i.e. phase postperiastron; False by default).

    Arguments
        phase : float
            Orbital phase to be computed.
        binary_system : binaries.System
            The binary system to be considered.
        initdate : float
            Minimum MJD epoch to compute.
        enddate : float
            Maximum MJD epoch to compute.
        postperiastron : bool [default: Fasle]
            In case the given phase is a post-periastron phase.
        superorbital : bool [default: False]
            In case the superorbital period (instead of the normal one) need to be computed.

    Returns
        times : list
            List of all MJD dates at which the system is the given orbital phase
            within the given temporal period.
    """
    t0, p, per = get_system_parameters(binary_system, superorbital=superorbital)
    # Determine a date at which the binary is in that phase:
    if postperiastron == True:
        phase = phase + per
    dayzero = t0 + phase*p
    # Just to avoid possible problems... (I should to check this and improve...)
    if dayzero > inidate:
        dayzero = dayzero - (np.ceil(dayzero/inidate))*p
    n_min = int(np.ceil( (inidate - dayzero)/p ))
    n_max = int(np.floor( (enddate - dayzero)/p ))
    return [ dayzero + n*p for n in range(n_min, n_max + 1) ]


if __name__ == '__main__':
    """If it is interactive mode, asks recursively a date and returns the orbital phase
    of the system on that date.
    NOTE: Have not touch the code since the very first version.
    """

    months = { 'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
    'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12  }
    ask_for_system = [True, True]
    binaries = binaries.Binaries()
    while ask_for_system[0] == True:
        # Ask for the binary system and take its parameters
        print('The available binary systems are:')
        binaries.list_binaries()
        binary = raw_input('Choose a binary system:\n')
        binary = binary.strip()
        try:
            T0, P, per = get_system_parameters(binary)
        except:
            answer = raw_input('Bad name. You want to continue? (y/n)\n')
            if answer == 'y' or answer == 'yes':
                continue
            else:
                break
        # date initially to zero (year, month, day, hour, minute, seccond)
        date = [0, 0, 0, 0, 0, 0]
        while ask_for_system[1] == True:
            thedate = ''
            try:
                # Ask for the date
                thedate = raw_input('Date in format: YYYY-Mmm-DD HH:MM:SS (HH:MM:SS optional)\n')

                # Get the year, the month and the day
                date[0] = int(thedate[0:4])
                date[1] = months[thedate[5:8]]
                date[2] = int(thedate[9:11])
                # Check if there are hours, minutes or secconds to include them
                length = len(thedate)
                # If there is only hours
                if length > 12:# and length < 15:
                    date[3] = int(thedate[12:14])
                    # If there is only hours and minutes
                    if length > 15:# and length < 18:
                        #date[3] = int(thedate[12:14])
                        date[4] = int(thedate[15:17])
                        # If there is hours, minutes and secconds
                        if length > 18:#: and length <= 20:
                            #date[3] = int(thedate[12:14])
                            #date[4] = int(thedate[15:17])
                            date[5] = int(thedate[18:20])
                dtdate = dt.datetime(date[0], date[1], date[2], date[3], date[4], date[5])
                phase = (mjd.date2mjd(dtdate) - T0)/P % 1
                if per != 0.0:
                    phase_postperiastron = (phase - per + 1) % 1
                    print('Orbital phase: ', phase)
                    print('Post-periastro phase: %.3f\n' % phase_postperiastron)
                else:
                    print('Orbital phase: %.3f\n' % phase)
            # If there is an exception. Because of a bad date or the user want to go out
            except:
                answer = raw_input('Bad date. Want continue? (y/n/o)\ny = yes, n = no, o = choose other binary\n')
                if answer == 'y' or answer == 'yes':
                    pass
                elif answer == 'o' or answer == 'other':
                    break
                else:
                    ask_for_system = [False, False]
                    break
                    #answer = raw_input('Want to choose other system? (y/n)\n')
                    #if answer == 'y' or answer == 'yes':
                    #    break
                    #else:
                    #    ask_for_system[0] = False
                    #    break


