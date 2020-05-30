"""
This module defines the class that deal with parameters.

It is unlikely users will need to work with these classes directly, unless they
define their own objects.

"""

from astropy import units as u

__author__ = "Benito Marcote"
__copyright__ = "Copyright 2014, Benito Marcote"
__credits__ = ["Benito Marcote"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Benito Marcote"
__email__ = "bmarcote@am.ub.es"
__status__ = "Development"


class Parameter(object):
    """
    Wraps individual parameters.

    This class represents a source's parameter (in a somewhat broad sense).  It
    acts as both a descriptor that can be assigned to a class attribute to
    describe the parameters accepted by an individual source.

    Parameter instances never store the actual value of the parameter
    directly.  Rather, each instance of a model stores its own parameters
    as either hidden attributes

    Parameters
    ----------
    name : str
        parameter name
    value : float or string
        value to use for this parameter
    description : str
        parameter description
    reference :str
        reference of the parameter value
    getter : callable
        a function that wraps the raw (internal) value of the parameter
        when returning the value through the parameter proxy (eg. a
        parameter may be stored internally as radians but returned to the
        user as degrees)
    setter : callable
        a function that wraps any values assigned to this parameter; should
        be the inverse of getter
    error : callable
        error in the parameter value
        Can be a string, a float or a 1-array with floats indicating the upper and lower errors.
    """

    def __init__(self, name, description='', value=None, has_error=False, error=None,
                reference=''):
        """Initialize a Parameter.
        If has_error is True, then the error attribute is ignored.
        """
        super(Parameter, self).__init__()
        self._name = name
        self.__doc__ = description.strip()
        self._value = value
        self._attr = '_'+name
        self._has_error = has_error
        self._error = error
        self._reference = reference

    def __getitem__(self):
        return self._value

    def __setitem__(self, value):
        self._value = value

    def __repr__(self):
        # if the value is given a a astropy (unit) Quantity with units...
        if type(self.value) == u.Quantity:
            return '{0!r}: {1!r} {2!r}'.format(self._name, self.value.value,
                    self.value.unit.name)
        else:
            return '{0!r}: {1!r}'.format(self._name, self.value)

    def __call__(self):
        return self._value

    @property
    def name(self):
        """Return the name of the parameter"""
        return self._name

    @property
    def description(self):
        return self.__doc__

    @description.setter
    def description(self, desc, append=True):
        if append:
            self._description = self._description + desc
        else:
            self._description = desc

    @property
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, reference, append=False):
        if append:
            self._reference = self._reference + reference
        else:
            self._reference = reference

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        self._value = val

    @property
    def error(self):
        if not self._has_error:
            raise AttributeError("This parameter has not errors")
        elif self._error is None:
            raise AttributeError("The error has not been setted.")
        else:
            return self._error

    @error.setter
    def error(self, val):
        if not self._has_error:
            raise AttributeError("This parameter has not errors")
        else:
            self._error = val

    def __add__(self, other):
        if type(other)==type(self):
            return self.value+other.value
        else:
            return self.value+other

    def __sub__(self, other):
        if type(other)==type(self):
            return self.value-other.value
        else:
            return self.value-other

    def __rsub__(self, other):
        if type(other)==type(self):
            return other.value-self.value
        else:
            return other-self.value

    def __mul__(self, other):
        if type(other)==type(self):
            return self.value*other.value
        else:
            return self.value*other

    def __div__(self, other):
        if type(other)==type(self):
            return self.value/other.value
        else:
            return self.value/other

    def __pow__(self, other):
        if type(other)==type(self):
            return self.value**other.value
        else:
            return self.value**other

    def __eq__(self, other):
        if type(other)==type(self):
            return self.value == other.value
        else:
            return self.value == other

    def __ne__(self, other):
        if type(other)==type(self):
            return not (self.value == other.value)
        else:
            return not (self.value == other)

    def __lt__(self, other):
        if type(other)==type(self):
            return self.value < other.value
        else:
            return self.value < other

    def __gt__(self, other):
        if type(other)==type(self):
            return self.value > other.value
        else:
            return self.value > other

    def __le__(self, other):
        if type(other)==type(self):
            return self.value <= other.value
        else:
            return self.value <= other

    def __ge__(self, other):
        if type(other)==type(self):
            return self.value >= other.value
        else:
            return self.value >= other


