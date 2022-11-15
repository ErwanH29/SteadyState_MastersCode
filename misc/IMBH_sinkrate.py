import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from amuse.ext.galactic_potentials import Plummer_profile, MWpotentialBovy2015, PowerLawCutoff_profile
from amuse.units import constants, units
import matplotlib.gridspec as gridspec

"""
Script to calculate the IMBH infall rate with many assumptions.
The data is taken from Ghez et al. 1999, 1997, Zhang and Fall 1999, Vudragovic et al. 2022,
Valenti 2016, POrtegeis Zwarte et al. 2006, Chandar et al. 2018

Assumptions made:
  -   ~10% of clusters form IMBH [Taken from Portegeis Zwart et al. 2006 and Hamilton and Miller 2002]
  -   Kroupa mass function of stars (Salpeter has a higher avg. mass and thus less stars per cluster and thus more clusters and more IMBH)
  -   Number of enclosed black hole assumes <m> ~ 0.4 though this may be larger for particles in the core
  -   Take the dynamical friction to sink all the way into the center, not to reach some distance r
"""

def tickers(ax):
    """
    Function to give outlay for the axis
    """

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.tick_params(axis="y", which = 'both', direction="in")
    ax.tick_params(axis="x", which = 'both', direction="in")    
    return ax

def IMBH_pop(MW_mass, clust_eff, avg_stellar):
  """
  Function to extract the ~ stellar population of the galactic core
  """

  no_stars = MW_mass/avg_stellar
  MW_clustS = no_stars * cluster_eff #number of stars from clusters
  no_clust = MW_clustS / no_clusterS
  no_IMBH = clust_eff * no_clust
  return no_IMBH

def df_timescale(dist, dispvel, IMBH_mass):
    return 1.9 * 10 **9 * ((dist)/5000)**2 * (dispvel)/200 * 10**8/(IMBH_mass)

#Units in MSun
MW_bulge = 2*10**10                           #Taken from Bovy 2015
avg_cluster = 10**5                           #Taken from TB2008
avg_stellarL = [0.4, 0.6]                     #Average of Kroupa and Salpeter mass DF
clust_effL = [0.05, 0.1]                      #Taken from Zhang and Fall 1999
cluster_eff = 0.25                            #Taken from Chandar et al. 2018
no_clusterS = avg_cluster/avg_stellarL[0]
IMBH_pop(MW_bulge, 0.2, 0.4)

plummer_prof = PowerLawCutoff_profile(2.22638e8|units.MSun/units.kpc**3, 1.|units.kpc, 1.8, 1.9|units.kpc)
swag = plummer_prof.enclosed_mass(8 | units.parsec)/10**6
distances = np.linspace(1,100,1000)
enc_mass = [4*10**6+(i/6.5)**(1.75)*(plummer_prof.enclosed_mass(i * 1 | units.parsec)).value_in(units.MSun) for i in distances] #Function based on Ghez et al. 1998
close_dist = np.linspace(0.1,10)
close_mass = [4*10**6 * i for i in close_dist]

plt.plot(distances, enc_mass)
plt.xscale('log')
plt.yscale('log')
plt.show()

no_IMBH = IMBH_pop(MW_bulge, clust_effL[0], avg_stellarL[0])
df_time = df_timescale(100, 50, 1000)
infall_rate = df_time / no_IMBH

df_timescale_list = [df_timescale(i, plummer_prof.circular_velocity(i | units.parsec).value_in(units.kms), 1000) for i in distances]
no_IMBH_list = [IMBH_pop(i, clust_effL[0], avg_stellarL[0]) for i in enc_mass]
IMBH_infall_rate = [i/(10**9) for i, j in zip (df_timescale_list, no_IMBH_list)]

df_timescale_list_eff = [df_timescale(i, plummer_prof.circular_velocity(i | units.parsec).value_in(units.kms), 1000) for i in distances]
no_IMBH_list_eff = [IMBH_pop(i, clust_effL[1], avg_stellarL[0]) for i in enc_mass]
IMBH_infall_rate_eff = [i/(10**9) for i, j in zip (df_timescale_list, no_IMBH_list)]

df_timescale_list_mass = [df_timescale(i, plummer_prof.circular_velocity(i | units.parsec).value_in(units.kms), 1000) for i in distances]
no_IMBH_list_mass = [IMBH_pop(i, clust_effL[0], avg_stellarL[1]) for i in enc_mass]
IMBH_infall_rate_mass = [i/(10**9) for i, j in zip (df_timescale_list, no_IMBH_list)]

idx = np.argwhere(np.asarray(IMBH_infall_rate) < 13.61)
distances = np.asarray(distances)
no_IMBH_list_cut = np.asarray(no_IMBH_list_mass)
IMBH_infall_rate = np.asarray(IMBH_infall_rate)

no_IMBH_list_Max = np.asarray([IMBH_pop(i, 3*clust_effL[1], avg_stellarL[0]) for i in enc_mass])
print('Max. Bound on infall rate: ', IMBH_infall_rate[idx][-1]/no_IMBH_list_cut[idx][-1] * 10**3)
print('Max. Bound on infall rate: ', IMBH_infall_rate[idx][-1]/no_IMBH_list_Max[idx][-1] * 10**3)

infall_rate = []
for i in range(len(idx)):
  val = IMBH_infall_rate[idx][i]/no_IMBH_list_cut[idx][i] * 10**3
  infall_rate.append(val)

infall_rate = np.asarray(infall_rate)
print(np.mean(infall_rate))
plt.plot(distances[idx], infall_rate)
plt.show()

plt.figure(figsize=(12, 4))
gs = gridspec.GridSpec(1, 2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax2.axhline(13.61, color = 'black', linestyle = '--')
ax2.text(1.2,17, 'Age of the MW')
ax2.plot(distances, IMBH_infall_rate, color = 'black')
ax1.plot(distances, no_IMBH_list, color = 'black',  label =r'$\langle m \rangle = 0.4, \epsilon_{IMBH} = 0.05$')
ax1.plot(distances, no_IMBH_list_eff, color = 'blue',  label =r'$\langle m \rangle = 0.4, \epsilon_{IMBH} = 0.10$')
ax1.plot(distances, no_IMBH_list_mass, color = 'red',  label =r'$\langle m \rangle = 0.6, \epsilon_{IMBH} = 0.05$')
ax1.set_ylim(0.1,5000)
for ax_ in [ax1, ax2]:
  ax_.set_xlim(1,100)
  ax_.set_xlabel('Distance to Core [pc]')
  tickers(ax_)
  ax_.set_yscale('log')
  ax_.set_xscale('log')

ax2.set_title('Distance to Core vs. IMBH Infall Rate')
ax1.set_title('Distance to Core vs. Enclosed IMBH Population')
ax2.set_ylabel(r'Dynamical Friction Time [Gyr]')
ax1.set_ylabel(r'Number of Enclosed IMBH')

ax1.legend(loc = 'upper left')
plt.savefig('test.pdf', dpi=300, bbox_inches='tight')
plt.clf()
