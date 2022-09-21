from amuse.ext.galactic_potentials import MWpotentialBovy2015
from parti_initialiser import *
import numpy as np
import matplotlib.pyplot as plt
import glob, os
import pickle as pkl

def potential_plotters():
    """
    Function to plot the potentials present in the simulation
    """
    
    SMBH_code = MW_SMBH()
    MWG_code  = MWpotentialBovy2015()

    dist_range = np.linspace(-0.1, 0.1, 10000) | units.parsec
    test_mass  = 100 | units.MSun
    r = 0.01 | units.parsec

    SMBH_potential = [(-test_mass*SMBH_code.get_potential_at_point(0, (i+r), (i+r), r*0.01**2)).value_in(units.J) for i in dist_range]
    MWG_potential  = [(-test_mass*MWG_code.get_potential_at_point(0 | units.kpc, (i+r), (i+r), r*0.01**2)).value_in(units.J) for i in dist_range]
    GC_potential   = [(-test_mass*GC_code.get_potential_at_point(0, i, i, r*0.01**2)).value_in(units.J) for i in dist_range]
    cum_potential  = [i+j+z for i,j,z in zip(GC_potential, SMBH_potential, MWG_potential)]

    plt.title('Potential Wells Used in the Simulation')
    plt.plot(dist_range.value_in(units.parsec), SMBH_potential, color = 'black', linestyle = '--', label = r'SMBH (m = $4\times10^6 M_{\odot}$)')
    plt.plot(dist_range.value_in(units.parsec), MWG_potential, color = 'black', linestyle = ':', label = 'Milky Way [Bovy (2015)]')
    plt.plot(dist_range.value_in(units.parsec), GC_potential, color = 'black', linestyle = '-.', label = 'Globular Cluster')
    plt.plot(dist_range.value_in(units.parsec), cum_potential, color = 'C0', label = 'Cumulative Potential')
    plt.xlabel(r'$R_c$ [parsec]')
    plt.ylabel(r'Binding Energy [Joules]')
    plt.legend()
    plt.xlim(-0.05, 0.05)
    plt.yscale('log')
    plt.savefig('figures/potentials_used.pdf', dpi = 300)

def find_nearest(array, value):
    """
    Function which finds the index of an element in an array
    corresponding to the nearest value of some input.
    
    Inputs:
    array:  Array for which we want to find the index of
    value:  Value for which we want to compare the array elements with
    output: Index where the value ~ array element
    """

    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index

def dynamical_fric(pos_distr, vel_distr, mass_distr, eta, tend):
    """
    Function which utilises eqn. 9 of Petts et al. 2012 to derive
    the dynamical friction [TO FIX]
    
    Inputs:
    pos_distr:  Array of some sample GC to compare positions (and enclosed mass) with
    vel_distr:  Array with some velocity distribution to find f(v<vi)
    mass_distr: Array of some sample GC to find enclosed mass
    eta:        Time-step used in the simulation
    tend:       The final time of the simulation
    output:     Dynamical friction in units of velocity
    """
    
    temp_pos   = np.linspace(0,0.4,1000) | units.parsec
    temp_mass  = 50 | units.MSun
    temp_mass2 = 200 | units.MSun
    temp_vel   = 5  | units.AU/units.yr
    temp_vel2  = 20  | units.AU/units.yr

    systv = (4/3)*np.pi*(1 | units.parsec)**3       # The volume for which the BHs are in (line 44 interface.py)
    systm = 10**7 | units.MSun                      # The mass of the GC which the BHs are in (line 44 interface.py)

    sample1 = [ ]
    sample2 = [ ]
    sample3 = [ ]

    for i in range(len(temp_pos)):
        index = find_nearest(pos_distr, temp_pos[i].value_in(units.AU))
        enc_mass = mass_distr[index] | units.MSun
        index2 = find_nearest(vel_distr, temp_vel.value_in(units.AU/units.yr))
        frac_vel = index2/len(vel_distr)

        value1 = (-16*np.pi**2*(constants.G)**2*(systm)*((systm)+temp_mass)*(systm/systv) \
                  *np.log((temp_mass/enc_mass)**2+1) *frac_vel*(temp_vel)**-2 * eta * tend).value_in(units.AU/units.yr)
        sample1.append(value1)
        

        value3 = (-16*np.pi**2*(constants.G)**2*(systm)*((systm)+temp_mass2)*(systm/systv) \
                  *np.log((temp_mass2/enc_mass)**2+1)*frac_vel*(temp_vel)**-2 * eta * tend).value_in(units.AU/units.yr)
        sample3.append(value3)
        
        index2 = find_nearest(vel_distr, temp_vel2.value_in(units.AU/units.yr))
        frac_vel = index2/len(vel_distr)
        value2 = (-16*np.pi**2*(constants.G)**2*(systm)*((systm)+temp_mass)*(systm/systv) \
                    *frac_vel*(temp_vel2)**-2 * eta * tend).value_in(units.AU/units.yr)
        sample2.append(value2)
        
    plt.plot(temp_pos.value_in(units.parsec), sample1, color = 'black',
                label = r'$M = 50M_{\odot}$, $v = 5.0$ AU yr$^{-1}$')
    plt.plot(temp_pos.value_in(units.parsec), sample3, color = 'red',
                ls = '--', label = r'$M = 200M_{\odot}$, $v = 5.0$ AU yr$^{-1}$')                
    plt.plot(temp_pos.value_in(units.parsec), sample2, color = 'blue', 
                ls = '--', label = r'$M = 50M_{\odot}$, $v = 20$ AU yr$^{-1}$')
                
    plt.title('Dynamical Friction vs.\nDistance from Cluster Core')
    plt.xlabel('Distance from Core [pc]')
    plt.ylabel('Dynamical Friction [AU/yr]')
    plt.legend()
    plt.show()
    #plt.savefig('figures/dynamicalfriction.pdf', dpi=300)
    return 

def velocityList():
    """
    Function to plot the Maxwellian distribution used
    """

    sigmaV = 6 # in kms
    vrange = np.linspace(0, 500, 10000) # in kms
    ProbFunc = [np.sqrt(2/np.pi)*(i**2/sigmaV**3)*np.exp(-i**2/(2*sigmaV**2)) for i in vrange]
    ProbFunc = ProbFunc/max(ProbFunc)
    CumSum = np.cumsum(ProbFunc)/200
    plt.plot(vrange, ProbFunc, color = 'red', label = 'PDF')
    plt.plot(vrange, CumSum, color = 'blue', label = 'CDF')
    plt.xlabel(r'Speed [km s$^{-1}$]')
    plt.ylabel(r'Normalised Population Fraction')
    plt.title('Velocity Distribution of Stars \n'
              'in Globular Clusters')
    plt.xlim(0,100)
    plt.legend()
    plt.savefig('figures/BasicMB.pdf', dpi = 300)

filename = glob.glob('data/preliminary_calcs/*')
with open(os.path.join(max(filename, key=os.path.getctime)), 'rb') as input_file:
    temp_data = pkl.load(input_file)
pos_values = temp_data.iloc[0][0]
mas_values = temp_data.iloc[0][1]
vel_values = temp_data.iloc[0][2]

temp_vel = dynamical_fric(pos_values, vel_values, mas_values, 10**-3, 100 | units.yr)