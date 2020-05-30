import numpy as np
from scipy.optimize import fsolve
from scipy import interpolate
from astropy import units as u

"""Let be a two-star system (A and B) with a wind collision region:

         CD
           \
A  <--Ra-- | <--Rb--  B
           /
 -------------D------

This program calculates the curve of the WCR (the contact discontinuity, CD)
where
   R      -- the distance from the star A to the WCR for an angle theta.
   Ra     -- the distance from the star A to the center of the WCR (in the line A-B).
   Rb     -- the distance from the star B to the center of the WCR (in the line A-B).
   theta  -- the angle from AB (theta=0) to the top following CD.
   theta1 -- the angle from B to CD, in an equivalent way to theta.
   D      -- the distance between A and B.
   eta = (Rb/Ra)**2  =  (Mdot_b * v_b) / (Mdot_a * v_a)

Using references Canto et al. (1996), Antokhin et al. (2004).
"""






class WCR(object):
    """Defines a Wind Collision Region (WCR)
    """
    def __init__(self, distance, eta):
        """Inputs:
            - distance : float or array-like
                Distance between the two stars.
            - eta : float
                The parameter eta = (Rb/Ra)**2  =  (Mdot_b * v_b) / (Mdot_a * v_a)
        """
        self._distance = distance
        self._eta = eta

    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self, new_distance):
        self._distance = new_distance

    @property
    def eta(self):
        return self._eta

    @eta.setter
    def eta(self, new_eta):
        self._eta = new_eta

    @property
    def Rb(self):
        """From Canto et al (1996):  Ra = sqrt(eta)*distance / (1 + sqrt(eta))
        And as D = Ra + Rb => Ra = sqrt(eta)*Rb
        """
        return self.distance/(1 + np.sqrt(self.eta))

    @property
    def Ra(self):
        """From Canto et al (1996):  Ra = sqrt(eta)*distance / (1 + sqrt(eta))
        """
        return np.sqrt(self.eta)*self.distance/(1 + np.sqrt(self.eta))

    @staticmethod
    def get_Ra_from_Rb(eta, Rb):
        """Determines Ra given eta and Rb
        """
        return Rb*np.sqrt(eta)

    @staticmethod
    def get_Rb_from_Ra(eta, Ra):
        """Determines Rb given eta and Ra
        """
        return Ra/np.sqrt(eta)

    @staticmethod
    def get_distance_from_Ra(eta, Ra):
        """Determines the distance between the two components given eta and Rb
        """
        return Ra*(1 + 1/np.sqrt(eta))

    @staticmethod
    def get_distance_from_Rb(eta, Rb):
        """Determines the distance between the two components given eta and Rb
        """
        return Rb*(1 + np.sqrt(eta))

    @property
    def theta_max(self):
        """Determines the maximum value for theta for a given CD with the parameter eta.
        """
        return fsolve(lambda t: t - np.tan(t) - np.pi/(1 - self.eta), np.pi/4)*u.rad

    @staticmethod
    def theta_max_from_eta(eta):
        """Determines the maximum value for theta for a given CD with the parameter eta.
        """
        return fsolve(lambda t: t - np.tan(t) - np.pi/(1 - eta), np.pi/4)*u.rad

    def get_theta1(self, theta):
        """Determines theta1 given eta and theta.
        If theta does not have units, radians are assumed.
        """
        if isinstance(theta, u.Quantity):
            f = lambda theta1 : np.nan_to_num(theta1/np.tan(theta1), nan=1.0) - 1 - \
                self.eta*(np.nan_to_num(theta.to(u.rad).value/np.tan(theta), nan=1.0) - 1)
            # return (fsolve(f, (np.nan_to_num(theta/np.abs(theta))).value)*u.rad).to(theta.unit)
            return (fsolve(f,  theta.to(u.rad).value)*u.rad).to(theta.unit)
        else:
            f = lambda theta1 : np.nan_to_num(theta1/np.tan(theta1), nan=1.0) - 1 - \
                self.eta*(np.nan_to_num(theta/np.tan(theta), nan=1.0) - 1)
            # return fsolve(f, np.zeros_like(theta1)+0.001)
            return fsolve(f, theta)

    def get_theta(self, theta1):
        """Determines theta given eta and theta
        If theta1 does not have units, radians are assumed.
        """
        if isinstance(theta1, u.Quantity):
            f = lambda theta : np.nan_to_num(theta1.to(u.rad).value/np.tan(theta1), nan=1.0) - 1 \
            - self.eta*(np.nan_to_num(theta/np.tan(theta), nan=1.0) - 1)
            # just to avoid to start looking for solutions in theta=0, as those angles do not converge
            # return (fsolve(f, np.zeros_like(theta1)+0.001*u.rad)*u.rad).to(theta1.unit)
            # return (fsolve(f, (theta1/np.abs(theta1)).value)*u.rad).to(theta1.unit)
            return (fsolve(f,  theta1.to(u.rad).value)*u.rad).to(theta.unit)
        else:
            f = lambda theta : np.nan_to_num(theta1/np.tan(theta1), nan=1.0) - 1 - \
                self.eta*(np.nan_to_num(theta/np.tan(theta), nan=1.0) - 1)
            return fsolve(f, theta1)

    def get_radius_from_theta_theta1(self, theta, theta1):
        """Determines the distance from A to CD for given angles theta and theta1
        Inputs
            theta  : float or astropy.Quantity
            theta1 : float or astropy.Quantity
            distance : float or astropy.Quantity
                distance between A and B.
        Returns
            radius : float or astropy.Quantity
                distance from A to the CD at the specified angles.
        """
        return self.distance*np.sin(theta1)/np.sin(theta + theta1)


    def get_radius_from_theta(self, theta):
        """Determines the distance from A to CD for a given angle theta
        assumis
            theta  : float or astropy.Quantity
        Returns:
            radius : float or astropy.Quantity
                distance from A to the CD at the specified angles.
        """
        theta1 = self.get_theta1(theta)
        radii = self.distance*np.nan_to_num(np.sin(theta1)/np.sin(theta + theta1))
        radii[np.where(theta == 0.0)] = self.Ra
        return radii

    def get_radius_from_theta1(self, theta1):
        """Determines the distance from A to CD for a given angle theta
        assumis
            theta1  : float or astropy.Quantity
        Returns:
            radius : float or astropy.Quantity
                distance from A to the CD at the specified angles.
        """
        theta = self.get_theta(theta1)
        radii = self.distance*np.nan_to_num(np.sin(theta1)/np.sin(theta + theta1))
        radii[np.where(theta == 0.0)] = self.Ra
        return radii

    def get_thetas(self, tolerance=15*u.deg):
        return np.linspace(-self.theta_max+tolerance, self.theta_max-tolerance, 100)[:,0]

    def get_radius(self, theta_tolerance=15*u.deg):
        return self.get_radius_from_theta(self.get_thetas(theta_tolerance))

    def get_x_y(self, theta_tolerance=15*u.deg):
        rs = self.get_radius(theta_tolerance)
        ts = self.get_thetas(theta_tolerance)
        return rs*np.cos(ts), rs*np.sin(ts)

    def get_x_y_from_thetas(self, thetas):
        rs = self.get_radius_from_theta(thetas)
        return rs*np.cos(thetas), rs*np.sin(thetas)


    def tangential_velocity_from_theta(self, theta, v_ratio):
        """Returns the tangential velocity, directed along the shell, of the flow
        Inputs:
            theta : float or astropy.Quantity
                The angle theta where to compute the velocity.
            v_ratio : float or astropy.Quantity (dimensionless)
                The ratio between wind velocities: v_a/v_b.
        """
        if isinstance(theta, u.Quantity):
            theta = theta.to(u.rad).value

        theta1 = self.get_theta1(theta)
        temp = np.sqrt((self.eta*(theta-np.sin(theta)*np.cos(theta)) + \
               (theta1-np.sin(theta1)*np.cos(theta1)))**2 + \
               (self.eta*np.sin(theta)**2-np.sin(theta1)**2)**2)
        return temp/(2*(self.eta*(1-np.cos(theta)) + v_ratio*(1-np.cos(theta1))))

    def mass_surface_density(self, theta, v_ratio, massloss_v_ratio):
        """Returns the mass surface density at the angle theta.
        Inputs:
            theta : float or astropy.Quantity
                The angle theta where to compute the velocity.
            v_ratio : float or astropy.Quantity (dimensionless)
                The ratio between wind velocities: v_a/v_b.
            massloss_v_ratio : float or astropy.Quantity
                Defined as Mdot_A / v_wind_A, where
                Mdot_A is the mass-loss rate of the star A.
                v_wind_A is the wind speed of the star A at the position of the shock.
        """
        sigma0 = massloss_v_ratio/(2*np.pi*self.eta*self.distance)
        if isinstance(theta, u.Quantity):
            theta = theta.to(u.rad).value

        theta1 = self.get_theta1(theta)
        temp = np.sin(theta+theta1)*(self.eta*(1-np.cos(theta))+v_ratio*(1-np.cos(theta1)))**2 \
                /(np.sin(theta)*np.sin(theta1))
        return sigma0*np.abs(temp)/np.sqrt(( self.eta*(theta-np.sin(theta)*np.cos(theta)) + \
            (theta1-np.sin(theta1)*np.cos(theta1)))**2 + \
            (self.eta*np.sin(theta)**2-np.sin(theta1)**2)**2 )


class WCR_rs(WCR):
    """Centers the stigmation point of the WCR in the (0,0) coordinates with a rotation of an angle alpha.
    The optional parameter y_is_declinationdetermines if the y axis refers to
    declination coordinates or not. In such case, in all transformations in the x axis,
    the declination will be taken into account.
    """
    def __init__(self, distance, eta, x0, y0, alpha, y_is_declination=False):
        self._x0 = x0
        self._y0 = y0
        self._alpha = alpha
        self.y_is_declination = y_is_declination
        super().__init__(distance, eta)

    @property
    def x0(self):
        return self._x0

    @x0.setter
    def x0(self, new_x0):
        self._x0 = new_x0

    @property
    def y0(self):
        return self._y0

    @y0.setter
    def y0(self, new_y0):
        self._y0 = new_y0

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, new_alpha):
        self._alpha = new_alpha

    def get_x_y(self):
        x, y = WCR.get_x_y(self)
        x2 = (x - self.Ra)*np.cos(self.alpha) - y*np.sin(self.alpha)
        y2 = (x - self.Ra)*np.sin(self.alpha) + y*np.cos(self.alpha)
        if self.y_is_declination:
            return x2/np.cos(self.y0) + self.x0, y2 + self.y0

        return x2 + self.x0, y2 + self.y0

    def get_x_y_from_thetas(self, thetas):
        x, y = WCR.get_x_y_from_thetas(self, thetas)
        x2 = (x - self.Ra)*np.cos(self.alpha) - y*np.sin(self.alpha)
        y2 = (x - self.Ra)*np.sin(self.alpha) + y*np.cos(self.alpha)
        if self.y_is_declination:
            return x2/np.cos(self.y0) + self.x0, y2 + self.y0

        return x2 + self.x0, y2 + self.y0

    def get_xyA(self):
        """Returns the x,y position of the star A.
        """
        x2 = -self.Ra*np.cos(self.alpha)
        y2 = -self.Ra*np.sin(self.alpha)
        if self.y_is_declination:
            return x2/np.cos(self.y0) + self.x0, y2 + self.y0

        return x2 + self.x0, y2 + self.y0

    def get_xyB(self):
        """Returns the x,y position of the star B.
        """
        x2 = self.Rb*np.cos(self.alpha)
        y2 = self.Rb*np.sin(self.alpha)
        if self.y_is_declination:
            return x2/np.cos(self.y0) + self.x0, y2 + self.y0

        return x2 + self.x0, y2 + self.y0


class WCR_2D(WCR_rs):
    """Allows to show the WCR as intensity in a 2D plane.
    Different functions can be assumed for how the intensity drops from the theta=0
    point. By defaults it assumes a gaussian profile.
    """

    @staticmethod
    def gaussian(theta, peak, sigma):
        return peak*np.exp(-0.5*(theta/sigma)**2)

    def griddata(self, xgrid, ygrid, profile, interpol='cubic', width=0.03, **kwargs):
        """Given a grid (xgrid, ygrid as created by np.meshgrid for two coordinate arrays)
        it returns the expected intensity for the given WCR.
        """
        x, y = self.get_x_y()
        z = profile(self.get_thetas(), **kwargs)
        # Creating values where intensity should be zero
        x2 = x - width*self.get_radius_from_theta(0.0)*np.cos(self.alpha)
        y2 = y - width*self.get_radius_from_theta(0.0)*np.sin(self.alpha)
        x3 = x + width*self.get_radius_from_theta(0.0)*np.cos(self.alpha)
        y3 = y + width*self.get_radius_from_theta(0.0)*np.sin(self.alpha)
        xn = np.append(x, x2)
        xn = np.append(xn, x3)
        yn = np.append(y, y2)
        yn = np.append(yn, y3)
        zn = np.append(z, np.zeros(len(x2) + len(x3)))
        zgrid = interpolate.griddata(np.array([xn, yn]).T, zn, np.array([xgrid, ygrid]).T,
                                     method=interpol, fill_value=0.0)
        return zgrid

    def get_z_from_grid(self, xgrid, ygrid, profile, **kwargs):
        """Given a grid (xgrid, ygrid as created by np.meshgrid for two coordinate arrays)
        it returns the expected intensity for the given WCR.
        """
        zgrid = np.zeros_like(xgrid)
        intensity = profile(self.get_thetas(), **kwargs)
        x,y = self.get_x_y()
        for xi,yi,zi in zip(x, y, intensity):
            i0, j0 = -1, -1
            for i in range(xgrid.shape[1]):
                if xgrid[0,i] <= xi:
                    i0 = i
                else:
                    break
            for j in range(ygrid.shape[0]):
                if ygrid[j,0] <= yi:
                    j0 = j
                else:
                    break
            zgrid[i0,j0] = zi
        return zgrid
        # when plotting, do zgrid.T as x,y axis are inverted from what I expected


class WCR_convolve(WCR_rs):
    """Convolves the WCR emission with the provided synthesized beam.
    """
    def __init__(self, distance, eta, x0, y0, alpha, beam_maj, beam_min, beam_pa):
        self._bmaj = beam_maj
        self._bmin = beam_min
        self._pa = beam_pa
        super().__init__(distance, eta, x0, y0, alpha)

    @property
    def synthesized_beam(self):
        return self._bmaj, self._bmin, self._bpa

    @synthesized_beam.setter
    def synthesized_beam(self, args):
        assert len(args) == 3
        self._bmaj, self._bmin, self._pa = args


    def get_z(self, xgrid, ygrid, peak, sigma, tolerance=0.001):
        """xgrid, ygrid must be a np.meshgrid of coordinates.
        """
        x, y = self.get_x_y()
        z = self.intensity_profile(x, peak, sigma)

        # DO THE OPPOSITE. RUN THROUGH x,y NAD GET THE CLOSEST X,Y point
        zgrid = np.zeros_like(xgrid)
        for i in range(zgrid.shape[0]):
            for j in range(zgrid.shape[1]):
                for k in range(len(x)):
                    if (np.abs(xgrid[i,j] - x[k]) < tolerance) and (np.abs(ygrid[i,j] - y[k]) < tolerance):
                        zgrid[i,j] = z[k]

        return zgrid


    def gaussian2d(self, xygrid, peak, sigma):
        x, y = xygrid
        return peak*np.exp( -((x-self.x0)**2 + (y-self.y0)**2)/(2*sigma**2) )


    # def gaussian2d_general(self.)


    def intensity_profile(self, coords, peak, sigma):
        """IMPORTANT. This function assumes a Gaussian profile for the WCR.
        The coords can be either the x values, the y ones, or the theta ones,
        as long as the WCR is centered (the stignation point is in the middle)
        of the array.
        sigma will have arbitrary units related to the length of coords
        """
        xs = np.arange(0.0, len(coords))
        xs -= (xs[-1]+xs[0])/2
        return peak*np.exp(-0.5*(xs/sigma)**2)




def hyperbola(x, eta, x0):
    b = np.arctan(max_theta(eta))*x0
    return b*np.sqrt((x/x0)**2 - 1)

def hyperbola_rs(x, eta, x0, y0, alpha):
    y = hyperbola(x, eta, x0)

