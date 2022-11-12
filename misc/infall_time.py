import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def IMBH_pop(MW_mass):
    no_stars = MW_mass/avg_stellar
    MW_clustS = no_stars * cluster_eff #number of stars from clusters
    no_clust = MW_clustS / no_clusterS
    no_IMBH = 0.1 * no_clust
    return no_IMBH

def df_timescale(dist, dispvel, IMBH_mass):
    return 1.9 * 10 **9 * (dist/5000)**2 * (dispvel)/200 * 10**8/(IMBH_mass)

#Units in MSun
MW_bulge = 2*10**10
avg_cluster = 10**5
avg_stellar = 0.4
cluster_eff = 0.25
no_clusterS = avg_cluster/avg_stellar

no_IMBH = IMBH_pop(MW_bulge)
df_time = df_timescale(100, 50, 1000)
infall_rate = df_time / no_IMBH
print('Infall Rate:', infall_rate)

distances = np.linspace(1,8,10000)
mass_distr = [4*10**6*i**0.75 for i in distances]
vel_dist  = [55/(np.log((i+0.01)))**0.05 for i in distances]

df_timescale_list = [df_timescale(i, j, 1000) for i, j in zip(distances, vel_dist)]
no_IMBH_list = [IMBH_pop(i) for i in mass_distr]
IMBH_infall_rate = [i/(10**6*j) for i, j in zip (df_timescale_list, no_IMBH_list)]

df_timescale_list_100 = [df_timescale(i, j, 100) for i, j in zip(distances, vel_dist)]
IMBH_infall_rate_100 = [i/(10**6*j) for i, j in zip (df_timescale_list_100, no_IMBH_list)]

fig, ax = plt.subplots()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.tick_params(axis="y", which = 'both', direction="in")
ax.tick_params(axis="x", which = 'both', direction="in")
plt.title('IMBH Infall Rate')
plt.plot(distances, IMBH_infall_rate, color = 'red',  label =r'$M_{IMBH} = 1000 M_{\odot}$')
plt.plot(distances, IMBH_infall_rate_100, color = 'blue', label =r'$M_{IMBH} = 100 M_{\odot}$')
plt.legend()
plt.show()