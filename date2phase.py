"""Convertions between dates and orbital phase for binary systems.
"""

import numpy as np
import datetime as dt
from astropy import units as u
import mjd
import .binaries as gb

def get_system_parameters(binary_system, superorbital=False):
    """ Returns the orbital parameters (T0, P and p.periastro) of the given binary systems.
    The valid binary systems are the binary systems available in Sources package."""
    # Copy the parameters from the imported file
    t0, p, phi = None, None, None
    try:
        obj = gb.Binaries().get_binary(binary_system)
        if superorbital:
            t0 = obj.t0_superorbital.value
            p = obj.period_superorbital.value.to(u.d).value
            phi = obj.periastron_phase.value
        else:
            t0 = obj.t0.value
            p = obj.period.value.to(u.d).value
            phi = 0.
        return t0, p, phi
    except KeyError:
        print('ERROR: System not found in the list or attributes T0, P or phi not found.')
        raise KeyError

    # Check if the object has the parameters



def date2phase(date, binary_system, superorbital=False):
    """ Return the orbital and post-periastron phase of the given binary systm for the input date.
    The date must be a datetime format."""
    MJD = mjd.date2mjd(date)
    pha, pha_pp = mjd2phase(MJD, binary_system, superorbital=superorbital)
    return pha, pha_pp

def mjd2phase(mjd, binary_system, superorbital=False):
    """ Return the orbital and post-periastron phase of the given binary systm for the input MJD.
    The MJD is the Modified Julian Date."""
    T0, P, per = get_system_parameters(binary_system, superorbital=superorbital)
    phase = (mjd - T0)/P % 1
    phase_postperiastron = phase
    if per != 0 and per != None:
        phase_postperiastron = (phase - per + 1) % 1
    return phase, phase_postperiastron

def phase2dates(phase, binary_system, inidate, enddate, postperiastron=False, superorbital=False):
    """ Find all the dates at which the binary system is in a determinate phase for a given range of dates.
    You can choose if that phase is or not a phase refered at the periastron (i.e. phase_postperiastron)
    by default that is False.
    Returns a list with these dates in MJD."""
    T0, P, per = get_system_parameters(binary_system, superorbital=superorbital)
    # Determine a date at which the binary is in that phase:
    if postperiastron == True:
        phase = phase + per
    dayzero = T0 + phase*P
    # Just to avoid possible problems... (I should to check this and improve...)
    if dayzero > inidate:
        dayzero = dayzero - (np.ceil(dayzero/inidate))*P
    n_min = int(np.ceil( (inidate - dayzero)/P ))
    n_max = int(np.floor( (enddate - dayzero)/P ))
    return [ dayzero + n*P for n in range(n_min, n_max+1) ]


if __name__ == '__main__':
    """If it is interactive mode, asks recursively a date and returns the orbital phase
    of the system on that date."""

    months = { 'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
    'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12  }
    ask_for_system = [True, True]
    binaries = gb.Binaries()
    while ask_for_system[0] == True:
        # Ask for the binary system and take its parameters
        print('Tha evailable binary systems are:')
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


