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