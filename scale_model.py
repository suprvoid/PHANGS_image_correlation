import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as s
from astropy.table import Table

class Scales_Model():


    def inputdata(self, radius, surf_den, v_disp, stell_den, v_circ, beta):

        self.radius = radius
        self.sigma = surf_den * s.M_sun.value/s.pc.value**2
        self.vdisp = v_disp * 1000
        self.rhost = stell_den * s.M_sun.value/s.pc.value**3
        self.V = v_circ * 1000
        self.beta = beta

        return
    

    def extractdata(self, filepath):

        table = Table.read(filepath)
        rad = np.array(table['r_gal'])
        surf_den = np.array(table['Sigma_mol'])
        v_disp = np.array(table['<vdisp_mol_pix_native>'])
        stell_den = np.array(table['rho_star_mp'])
        v_circ = np.array(table['V_circ_CO21_lgd'])
        beta = np.array(table['beta_CO21_lgd'])

        self.inputdata(rad, surf_den, v_disp, stell_den, v_circ, beta)

        return
    

    def scaleheight(self):

        valid_i = np.where(~np.isnan(self.rhost) & ~np.isnan(self.sigma) & ~np.isnan(self.vdisp))[0]

        quadratic = np.array([(4*np.pi*s.G.value*self.rhost[valid_i]), (np.sqrt(8*np.pi)*s.G.value*self.sigma[valid_i]), -self.vdisp[valid_i]**2])

        N = quadratic.shape[1]

        H = np.zeros(N)
        for i in range(N):

            roots = np.roots(quadratic[:,i])
            proot = roots[roots>0][0]
            H[i] = proot

        H /= s.pc.value

        self.H = H
        self.H_rad = self.radius[valid_i]

        return H
    

    def wave2dscale(self):

        valid_i = np.where(~np.isnan(self.vdisp) & ~np.isnan(self.sigma))[0]

        lamb2d = (2*self.vdisp[valid_i]**2)/(s.G.value*self.sigma[valid_i])
        lamb2d /= s.pc.value

        self.lamb2d = lamb2d
        self.lamb2d_rad = self.radius[valid_i]

        return lamb2d

    
    def jeanscale(self):
        
        self.scaleheight()

        valid_i = np.where(~np.isnan(self.sigma))[0][:self.H.size]

        self.rho = self.sigma[valid_i]/(np.sqrt(2*np.pi)*self.H*s.pc.value)

        jean = (self.vdisp[valid_i]*np.pi**(1/2))/(s.G.value*self.rho)**(1/2)
        jean /= s.pc.value

        self.jean = jean
        self.jean_rad = self.radius[valid_i]

        return jean


    def toomrescale(self):

        valid_i = np.where(~np.isnan(self.beta) & ~np.isnan(self.sigma))[0]

        self.kappa = np.sqrt(2*(1+self.beta[valid_i]))*self.V[valid_i]/(self.radius[valid_i]*1000*s.pc.value)

        toomre = 4*(np.pi**2)*s.G.value*self.sigma[valid_i]/(self.kappa**2)
        toomre /= s.pc.value

        self.toomre = toomre
        self.toomre_rad = self.radius[valid_i]

        return toomre
    

    def allscales(self):

        self.scaleheight()
        self.wave2dscale()
        self.jeanscale()
        self.toomrescale()

        return

    
    def plot_allscales(self, measurements, uncertainty):

        alp = 0.3
        siz = 5

        errors = np.log10(1+uncertainty/measurements[0])*np.log10(measurements[0])

        plt.figure()
        plt.style.use('bmh')
        plt.plot(self.H_rad, np.log10(self.H), marker='v', c='m', label='scale height', alpha=alp, markersize=siz)
        plt.plot(self.lamb2d_rad, np.log10(self.lamb2d), marker='v', c='g', label='most unstable 2D', alpha=alp, markersize=siz)
        plt.plot(self.jean_rad, np.log10(self.jean), marker='v', c='r', label='Jeans length', alpha=alp, markersize=siz)
        plt.plot(self.toomre_rad, np.log10(self.toomre), marker='v', c='b', label='Toomre length', alpha=alp, markersize=siz)
        plt.errorbar(measurements[1], np.log10(measurements[0]), xerr=None, yerr=errors, lw=1.6, fmt='X', c='k', label='measured scales', markersize=7)
        plt.ylabel('log \u03BB (pc)')
        plt.xlabel('Galactic radius (kpc)')
        plt.title(self.name)
        plt.legend(fontsize=8, loc='lower right')
        plt.show()
        plt.style.use('default')

        return


    def inst(self, name, filepath):

        self.name = name

        self.extractdata(filepath)

        self.allscales()

        return