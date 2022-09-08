
from amuse.lab import *
import matplotlib.pyplot as plt
import numpy as np

n_stars = 10**8

m_stars2 = new_powerlaw_mass_distribution(n_stars, 1|units.MSun, 
                                          100.0|units.MSun, alpha = -2.35)

xmin = np.log10(0.5)
xmax = np.log10(150)
bins = 10**(np.linspace(xmin, xmax, 51))
Nbin2, bin_edges2 = np.histogram(m_stars2.value_in(units.MSun), bins=bins)

y_stars2 = Nbin2 / (bin_edges2[1:] - bin_edges2[:-1])
x_stars  = (bin_edges2[1:] + bin_edges2[:-1]) / 2.0

for i in range(len(y_stars2)):
    y_stars2[i] = max(y_stars2[i], 1.e-10)

c2 = 0.005*((m_stars2.max().value_in(units.MSun)**(-2.35+1)) \
         	- (m_stars2.min().value_in(units.MSun)**(-2.35+1))) / (-2.35+1)

plt.title('Population Distribution for the Salpeter Model')
plt.fill_between(x= x_stars, 
                 y1= 9000 * len(y_stars2)/c2 * (x_stars**-2.35), 
                 where= x_stars > 10, color= "b", alpha= 0.2)
plt.plot(x_stars, 9000 * len(y_stars2)/c2 * (x_stars**-2.35), color = 'black', linestyle = '--')
plt.scatter(x_stars, y_stars2,s=20, color = 'C0', label = r'$\alpha = -2.35$')
plt.xlim(0.5, 150)
plt.ylim(1, 10**10)
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'Population [$N$]')
plt.xlabel(r'Mass [$M_\odot$]')
plt.legend()
plt.savefig('plummer_stats/plummer_dist.pdf', dpi = 300)
