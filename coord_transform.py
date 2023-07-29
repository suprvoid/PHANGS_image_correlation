import numpy as np
import astropy.units as u
import astropy.coordinates as c


class Galaxy:

    def setpar(self, ra_in, dec_in, incl_in, PA_in, dist_in):

        self.coords = c.SkyCoord(ra_in, dec_in, unit='deg')
        self.incl = incl_in*u.deg
        self.PA = PA_in*u.deg
        self.dist = dist_in

        return


def ang_coord(Galaxy, SkyCoord):
    
    phi = Galaxy.coords.position_angle(SkyCoord).degree-Galaxy.PA.value
    
    return phi*u.deg


def polar(Galaxy, SkyCoord):

    angle = c.angular_separation(Galaxy.coords.ra, Galaxy.coords.dec, SkyCoord.ra, SkyCoord.dec).value
    phi = ang_coord(Galaxy, SkyCoord)
    R = np.tan(angle)*Galaxy.dist*(np.cos(phi)**2+(np.sin(phi)/np.cos(Galaxy.incl))**2)**(1/2)

    return [R, phi]

def cartesian(Galaxy, polarvec):

    # R in kpc, phi in deg
    R, phi = polarvec[0], polarvec[1]
    L = (R/((np.cos(phi)**2+(np.sin(phi)/np.cos(Galaxy.incl))**2)**(1/2)*Galaxy.dist))*180/np.pi
    PA = phi+Galaxy.PA
    ra = Galaxy.coords.ra.value+L*np.sin(PA)
    dec = Galaxy.coords.dec.value+L*np.cos(PA)

    return [ra, dec]