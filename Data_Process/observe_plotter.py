from spatial_plotters import *
import numpy as np
import matplotlib.pyplot as plt

# Data is taken from 2022MNRAS.509.1587 [Walton et al. 2022] and provided by VizieR

data = np.loadtxt('data/observational_data/asu.tsv', dtype = str, delimiter='|')

rasc = []       #In deg
decl = []       #In deg
gal_name = []
gal_dist = []   #In Mpc
gal_type = []
avg_flux = []   #In mW/m2
peak_flux = []
gal_sep = []    #In arcsec
avg_lum = []    #In 1e-7 W
peak_lum = []   

for row_ in data[3:]:
    for i in range(np.shape(data[3:])[1]-1):
        row_[i].strip()
    
    if '1' in row_[9]:
        rasc.append(row_[0])
        decl.append(row_[1])
        gal_name.append(row_[2])
        gal_dist.append(row_[3])
        gal_type.append(row_[4])
        avg_flux.append(row_[6])
        peak_flux.append(row_[7])
        gal_sep.append(row_[8])
        avg_lum.append(row_[9])
        peak_lum.append(row_[10])

    elif '1' in row_[14]:
        rasc.append(row_[0])
        decl.append(row_[1])
        gal_name.append(row_[2])
        gal_dist.append(row_[3])
        gal_type.append(row_[4])
        avg_flux.append(0)
        peak_flux.append(row_[12])
        gal_sep.append(row_[13])
        avg_lum.append(0)
        peak_lum.append(row_[14])

    elif '1' in row_[18]:
        rasc.append(row_[0])
        decl.append(row_[1])
        gal_name.append(row_[2])
        gal_dist.append(row_[3])
        gal_type.append(row_[4])
        avg_flux.append(0)
        peak_flux.append(row_[16])
        gal_sep.append(row_[17])
        avg_lum.append(0)
        peak_lum.append(row_[18])
        
peak_lum = np.asarray([float(i) for i in peak_lum])
gal_dist = np.asarray([float(i) for i in gal_dist])[peak_lum > 10**40]
avg_flux = np.asarray([float(i) for i in avg_flux])[peak_lum > 10**40]
peak_flux = np.asarray([float(i) for i in peak_flux])[peak_lum > 10**40]
gal_sep = np.asarray([float(i) * 4.84814*10**-6 for i in gal_sep])[peak_lum > 10**40]  #Convert arcsec to rad
gal_sep = np.asarray([dist_*10**6 * np.arctan(ang_) for dist_, ang_ in zip(gal_dist, gal_sep)])
avg_lum = np.asarray([float(i) for i in avg_lum])[peak_lum > 10**40]
peak_lum = peak_lum[peak_lum > 10**40]
print(np.shape(gal_sep), np.shape(peak_lum))

plot_ini = plotter_setup()

fig, ax = plt.subplots()
colour_axes = ax.scatter(np.log10(gal_sep), np.log10(peak_lum), c = np.log10(avg_flux), edgecolors = 'black')
plot_ini.tickers(ax, 'plot')
ax.set_xlabel('Galactic Separation [pc]')
ax.set_ylabel(r'Peak Luminosity [erg s$^{-1}$]')
ax.set_xlim(0, 1.05*max(np.log10(gal_sep)))
plt.colorbar(colour_axes, ax = ax, label = r'$\log_{10}\langle f\rangle$ [mW m$^{-2}$]')
plt.show()