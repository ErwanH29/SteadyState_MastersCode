from amuse.ext.galactic_potentials import MWpotentialBovy2015
from initialiser import *
import numpy as np
import matplotlib.pyplot as plt

SMBH_code = MW_SMBH()
MWG_code  = MWpotentialBovy2015()
GC_code   = GC_pot()

dist_range = np.linspace(-0.1, 0.1, 10000) | units.parsec
test_mass  = 100 | units.MSun
r = 0.01 | units.parsec

SMBH_potential = [(-test_mass*SMBH_code.get_potential_at_point(0, (i+r), (i+r), r*0.01**2)).value_in(units.J) for i in dist_range]
MWG_potential  = [(-test_mass*MWG_code.get_potential_at_point(0 | units.kpc, (i+r), (i+r), r*0.01**2)).value_in(units.J) for i in dist_range]
GC_potential   = [(-test_mass*GC_code.get_potential_at_point(0, i, i, r*0.01**2)).value_in(units.J) for i in dist_range]
cum_potential = [i+j+z for i,j,z in zip(GC_potential, SMBH_potential, MWG_potential)]

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