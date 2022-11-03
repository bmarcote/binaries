#! /usr/bin/env python3
"""
Physical and orbital properties of the known gamma-ray binaries
and some X-ray binaries.

"""

import os as _os
import configparser as _configparser
from astropy import units as u
from astropy import coordinates
from astropy import time
from binaries.src.parameters import Parameter

__author__ = "Benito Marcote"
__copyright__ = "Copyright 2014, Benito Marcote"
__credits__ = ["Benito Marcote"]
__license__ = "GPL"
__version__ = "2.0.1"
__maintainer__ = "Benito Marcote"
__email__ = "bmarcote@am.ub.es"
__status__ = "Development"


class System(object):
    def __init__(self, name, coordinates, **kwargs):
        self._name = name
        self._coordinates = coordinates
        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def name(self):
        return self._name

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates):
        self._coordinates = coordinates

    def __repr__(self):
        return '<Sources.binaries.System: '+self.name+'>'


def read_sources_from_inifiles(path):
    _config = _configparser.ConfigParser()
    _binaries = {}
    # Import all the known objects
    _files = _os.listdir(path)
    for a_file in _files:
        if a_file[-4:] == '.ini':
            _config.read(path+'/'+a_file)
            #if is_python27:
            #    _config = _config.__dict__['_sections']
            #
            params = {}
            source_name = ''
            abrv_name = ''
            for a_section in _config.sections():
                if a_section == 'source_name':
                    source_name = _config[a_section]['name']
                    abrv_name = _config[a_section]['abrv']
                    continue

                name = a_section
                desc = _config[a_section].get('description', '')
                ref = _config[a_section].get('reference', '')
                val = _config[a_section].get('value')
                err = _config[a_section].get('error', None) # Change this
                has_error = _config[a_section].getboolean('has_error', fallback=True)
                try:
                    val = float(val)
                    err = [float(i) for i in err.split(',')]
                    #err = float(err)
                    if _config.has_option(a_section, 'unit') and (name != 'ra' or name != 'dec'):
                        val = val*eval(_config[a_section]['unit'])
                        if err != None:
                            # Can be two values: -err1, +err2 or just one +-err
                            for i in range(len(err)):
                                err[i] = err[i]*eval(_config[a_section]['unit'])

                    # If there is only one value for the error, then ocnvert to a quantity.
                    if len(err) == 1:
                        err = err[0]

                except (TypeError, ValueError):
                    # The values are strings
                    pass

                params[name] = Parameter(name, desc, val, has_error, err, ref)

            ra = params['ra']
            dec = params['dec']
            params.pop('ra')
            params.pop('dec')
            # ra and dec are coordinates
            ra_units = eval(_config['ra']['unit'])
            dec_units = eval(_config['dec']['unit'])
            distance = params['distance']
            if 'mu_alpha_cos_delta' in params:  # It has proper motions defined!
                pm_ra_cos_dec = params['mu_alpha_cos_delta']
                pm_dec = params['mu_delta']
                if 'ref_epoch' in params:
                    ref_epoch = time.Time(params['ref_epoch'].value, format='mjd')
                else:
                    print('WARNING: Epoch J2000.0 is assumed for {}.'.format(source_name))
                    ref_epoch = time.Time(2000.0, format='decimalyear')

                coord_val = coordinates.SkyCoord(ra=ra.value, dec=dec.value, unit=(ra_units,dec_units),
                                                 distance=distance.value, pm_ra_cosdec=pm_ra_cos_dec.value,
                                                 pm_dec=pm_dec.value, obstime=ref_epoch)
                coord_err = coordinates.SkyCoord(ra=ra.error, dec=dec.error, unit=(ra_units, dec_units),
                                                 distance=distance.error)
                                                 #pm_dec=pm_dec.value,
                                                 #pm_ra_cosdec=pm_ra_cos_dec.value, obstime=ref_epoch)
            else:
                coord_val = coordinates.SkyCoord(ra=ra.value, dec=dec.value, unit=(ra_units,dec_units),
                                                 distance=distance.value)
                coord_err = coordinates.SkyCoord(ra=ra.error, dec=dec.error, unit=(ra_units, dec_units),
                                                 distance=distance.value)

            ref = ra.reference.strip()
            if ref != dec.reference.strip():
                ref = ','.join([ref, dec.reference.strip()])
            coord = Parameter('Coordinates', ra.description+'\n'+dec.description, coord_val, True,
                coord_err, ref)
            _binaries[abrv_name] = System(source_name, coord, **params)
    return _binaries


class Binaries(object):
    """
    The initialized object contains all the known sources with its properties as attributes.
    """
    def __init__(self):
        path = _os.path.abspath(__file__)
        dir_path = _os.path.dirname(path)
        self._binaries = read_sources_from_inifiles(dir_path + '/data/')
        for key in self._binaries.keys():
            #print(key)
            setattr(self, key, self._binaries[key])

    def list_binaries(self):
        for key in self._binaries.keys():
            print(self._binaries[key].name)

    def return_name_binaries(self):
        names = []
        for key in self._binaries.keys():
            names.append(self._binaries[key].name)

        return names

    def get_binaries_as_dict(self):
        d = {}
        for key in self._binaries.keys():
            d[self._binaries[key].name] = self._binaries[key]

        return d

    def get_binary(self, name):
        for key in self._binaries:
            if self._binaries[key].name == name:
                return self._binaries[key]

        # No binary found
        raise KeyError

    def reload(self):
        """
        Reload all the binaries.
        Read again the existing init files to load all the known systems.
        The new object has all the sytems and their parameters updated.
        """
        # First remove all the current attributes
        for key in _binaries.keys():
            delattr(self, key)

        # Read again the init files and put the new attributes
        path = _os.path.abspath(__file__)
        dir_path = _os.path.dirname(path)
        self._binaries = read_sources_from_inifiles(dir_path + '/data/')
        for key in _binaries.keys():
            print(key)
            setattr(self, key, self._binaries[key])



