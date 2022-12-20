import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def MB_distr(vel):
    sigmaV = 150
    return np.sqrt(2/np.pi)*(vel**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

fig, ax = plt.subplots()
vel_list = np.linspace(0,700,1000)
vel_prob = [MB_distr(i) for i in vel_list]
vel_prob /= max(vel_prob)
ax.plot(vel_list, vel_prob, color = 'black')
plt.title('Velocity Distribution')
ax.set_xlabel(r'$|v|$ [km/s]')
ax.set_ylabel(r'$P(v)/P(v)_{\rm{max}}$')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.tick_params(axis="y", which = 'both', direction="in")
ax.tick_params(axis="x", which = 'both', direction="in")
plt.savefig('figures/velocity_distribution.pdf',  dpi=500, bbox_inches='tight')
plt.clf()

#Cluster mass distribution from Rodriguez et al. 2016a

mu, sigma = np.log(10**5.54), 0.52
s = np.random.lognormal(mu, sigma, 1000000)


count, bins, ignored = plt.hist(s, 100, align='mid', color = 'black', alpha = 0.3)
plt.clf()

fig, ax = plt.subplots()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
x = np.linspace(1.2*min(bins), 10*max(bins), 10000)

pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
       / (x * sigma * np.sqrt(2 * np.pi)))
x = [np.log10(4.5*i/10**5) for i in x]

pdf /= max(pdf)

ax.set_xlabel(r'$\log_{10}M_{\rm{GC}}$ [$\frac{M_{\rm{GC}}}{10^{5}M_\odot}$]')
ax.set_ylabel(r'$\rho/\rho_{\rm{max}}$')
ax.tick_params(axis="y", which = 'both', direction="in")
ax.tick_params(axis="x", which = 'both', direction="in")
ax.plot(x, pdf, linewidth=2, color='black')
ax.axis('tight')
plt.savefig('figures/globular_cluster_mass_distr.pdf',  dpi=500, bbox_inches='tight')

k = np.sort(bins)
print(len(k[5*k>10**6])/len(k))