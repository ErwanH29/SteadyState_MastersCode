import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def MB_distr(vel):
    sigmaV = 150
    return np.sqrt(2/np.pi)*(vel**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

fig, ax = plt.subplots()
vel_list = np.linspace(0,700,1000)
vel_prob = [MB_distr(i) for i in vel_list]
ax.plot(vel_list, vel_prob, color = 'black')
plt.title('Velocity Distribution')
ax.set_xlabel('Velocity [km/s]')
ax.set_ylabel('Probability')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.tick_params(axis="y", which = 'both', direction="in")
ax.tick_params(axis="x", which = 'both', direction="in")
plt.savefig('figures/velocity_distribution.pdf',  dpi=500, bbox_inches='tight')