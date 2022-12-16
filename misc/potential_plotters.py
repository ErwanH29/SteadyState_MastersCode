from amuse.ext.galactic_potentials import MWpotentialBovy2015
from amuse.lab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def potential_plotters():
    """
    Function to plot the potentials present in the simulation
    """
    
    MWG_code  = MWpotentialBovy2015()

    dist_range = np.linspace(-2, 2, 10000) | units.parsec
    test_mass  = 1000 | units.MSun
    r = 0.1 | units.parsec

    SMBH_potential = [(-test_mass*(-constants.G*(4*10**6 | units.MSun))/(abs(i))).value_in(units.J) for i in dist_range]
    MWG_potential  = [(-test_mass*MWG_code.get_potential_at_point(0 | units.kpc, (i+r), (i+r), r*0.01**2)).value_in(units.J) for i in dist_range]    
    cum_potential  = [i+j for i,j in zip(SMBH_potential, MWG_potential)]

    #MW domiantes potential at r ~ 0.002pc. The esc. velocity for this distance is:
    vesc = np.sqrt(2)*np.sqrt(abs(MWG_code.get_potential_at_point(0 | units.kpc, (r), (r), r))).in_(units.pc/units.yr)
    print('Escape vel.: ', vesc)
    #Thus the crossing time for my ejection condition is:
    tcross = ((2 | units.parsec)/vesc).in_(units.yr)
    print('Cross time', tcross, '. \nThis is: ', (tcross/1000).value_in(units.yr), ' timesteps.')

    fig, ax = plt.subplots()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.tick_params(axis="y", which = 'both', direction="in")
    ax.tick_params(axis="x", which = 'both', direction="in") 

    ax.set_title('Milky Way Core Potential Wells')
    ax.plot(dist_range.value_in(units.parsec), SMBH_potential, color = 'blue', label = r'SMBH (m = $4\times10^6 M_{\odot}$)')
    ax.plot(dist_range.value_in(units.parsec), MWG_potential, color = 'red', label = 'Milky Way [Bovy (2015)]')
    ax.plot(dist_range.value_in(units.parsec), cum_potential, color = 'black', label = 'Cumulative Potential')
    ax.set_xlabel(r'$r$ [pc]')
    ax.set_ylabel(r'$|P_E|$ [J]')
    ax.set_xlim(3e-4, 2)
    plt.legend()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.savefig('figures/potentials_used.pdf', dpi = 300, bbox_inches='tight')

def find_nearest(array, value):
    """
    Function which finds the index of an element in an array
    corresponding to the nearest value of some input.
    
    Inputs:
    array:  Array for which we want to find the index of
    value:  Value for which we want to compare the array elements with
    """

    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index

def velocityList():
    """
    Function to plot the Maxwellian distribution used.
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

potential_plotters()
