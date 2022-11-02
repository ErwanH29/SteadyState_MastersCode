from amuse.lab import *
from spatial_plotters import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

class bin_tert_systems(object):
    """
    Class which forms plots of the sustainable hierarchical and binaries 
    along with final time-step binary/hierarchical data.
    Sustainable condition: If same partner for iter > 5, or roughly 5000 years
    """

    def __init__(self):
        """
        Extracts the required data
        """
        
        np.seterr(divide='ignore')

        H0 = 73.04 #Taken from arXiv:2112.04510 in km/s/Mpc
        self.tH = (H0*(3.2408*10**-20))**-1

        self.IMBH_tracker = bulk_stat_extractor('data/Hermite/particle_trajectory/*', 'Y')
        self.no_data = len(self.IMBH_tracker)
    
        self.init_pop = []
        self.semi_SMBH_fin = []
        self.semi_SMBH_ini = []
        self.ecc_SMBH_fin = []
        self.ecc_SMBH_ini = []
        self.close_enc = []
        self.mass = []

        self.semi_bin_SMBH = []
        self.semi_bin_fin = []
        self.ecc_bin_SMBH = []
        self.ecc_bin_fin = []
        self.mass_bin_SMBH = []
        self.mass_bin = []
        
        for j in range(self.no_data):
            sim_data = self.IMBH_tracker[j]
            for parti_ in range(len(sim_data)):
                Nenc = 0

                if parti_ == 0:  # Don't care for SMBH
                    pass
                else:
                    if isinstance(sim_data.iloc[parti_][-1][0], np.uint64):           # Neglect removed particle
                        self.semi_SMBH_fin.append(sim_data.iloc[parti_][-2][7][0])    # Preserve SMBH data 
                        self.semi_SMBH_ini.append(sim_data.iloc[parti_][2][7][0]) 
                        self.ecc_SMBH_fin.append(1-sim_data.iloc[parti_][-2][8][0])
                        self.ecc_SMBH_ini.append(1-sim_data.iloc[parti_][2][8][0])
                        self.mass.append([sim_data.iloc[0][0][1], sim_data.iloc[parti_][0][1]])
                        
                        val = len(sim_data)
                        if val % 10 != 0:
                            val -=1
                        self.init_pop.append(val)

                        if sim_data.iloc[parti_][-2][8][1] < 1 and abs(sim_data.iloc[parti_][-2][7][1]) < 0.15 | units.parsec:    #Second Nearest binary and only look if eccentricity < 1
                            if sim_data.iloc[parti_][-2][6][1] == sim_data.iloc[0][0][0]:                    #Do not consider SMBH
                                pass
                            else:
                                print('Binary detected!')
                                self.semi_bin_SMBH.append(sim_data.iloc[parti_][-2][7][0])
                                self.semi_bin_fin.append(sim_data.iloc[parti_][-2][7][1])
                                self.ecc_bin_SMBH.append(1-sim_data.iloc[parti_][-2][8][0])
                                self.ecc_bin_fin.append(1-sim_data.iloc[parti_][-2][8][1])
                                self.mass_bin_SMBH.append([sim_data.iloc[0][0][1], sim_data.iloc[parti_][0][1]])
                                for part_ in range(len(sim_data)):
                                    if sim_data.iloc[parti_][-2][6][1] == sim_data.iloc[part_][0][0]:
                                        m2 = sim_data.iloc[part_][0][1]
                                self.mass_bin.append([sim_data.iloc[0][0][1], m2])

                        for k in range(np.shape(sim_data)[1]):
                            if sim_data.iloc[parti_][k][-1] < 5e-3:    # 5e-3 pc ~ 1000 AU and corresponds to ~5e42 J, a large KE when looking at their typical values
                                Nenc += 1/2
                        self.close_enc.append(Nenc)
    
    def gw_calc(self, semi, ecc, m1, m2):
        red_mass = (m1*m2)/(m1+m2)
        tot_mass = m1 + m2
        tgw = (5/256) * (constants.c)**5/(constants.G**3)*(semi**4*(1-ecc**2)**3.5)/(red_mass*tot_mass**2)
        return tgw

    def semi_ecc_SMBH_gw_plotter(self):

        plot_init = plotter_setup()
        tgw_SMBH_fin = []
        tgw_SMBH_ini = []
        ecc_evol = []
        
        for i in range(len(self.semi_SMBH_fin)):
            grav_ftimescale = self.gw_calc(self.semi_SMBH_fin[i], self.ecc_SMBH_fin[i], self.mass[i][0], self.mass[i][1]).value_in(units.s)
            grav_itimescale = self.gw_calc(self.semi_SMBH_ini[i], self.ecc_SMBH_ini[i], self.mass[i][0], self.mass[i][1]).value_in(units.s)
            tgw_SMBH_fin.append(grav_ftimescale)
            tgw_SMBH_ini.append(grav_itimescale)
        tgw_SMBH_fin = np.asarray(tgw_SMBH_fin)
        tgw_SMBH_ini = np.asarray(tgw_SMBH_ini)
        tgw_evol_SMBH = np.asarray([i/j for i, j in zip(tgw_SMBH_fin, tgw_SMBH_ini)])
        ecc_evol = [i/j for i, j in zip(self.ecc_SMBH_fin, self.ecc_SMBH_ini)]
        semi_evol = [i/j for i, j in zip(self.semi_SMBH_fin, self.semi_SMBH_ini)]

        tcomp = np.asarray([i/self.tH for i in tgw_SMBH_fin])
        semi_SMBH = np.asarray([i.value_in(units.parsec) for i in self.semi_SMBH_fin])
        ecc_SMBH  = np.asarray([np.log10(i) for i in np.asarray(self.ecc_SMBH_fin)])

        tcomp_crop = tcomp[tcomp<1e9]
        index_excess_crop = np.where(tcomp_crop > 1)[0]
        index_merger_crop = np.where(tcomp_crop <= 1)[0]
        tgw_SMBH_fin_crop = tgw_SMBH_fin[tcomp<1e9]
        tgw_SMBH_ini_crop = tgw_SMBH_ini[tcomp<1e9]
        tgw_SMBH_fin_ratio = tgw_SMBH_fin / self.tH
        tgw_SMBH_ini_ratio = tgw_SMBH_ini / self.tH
        
        exc_norm_tgw_fin_crop = tgw_SMBH_fin_crop[index_excess_crop] / max(tgw_SMBH_ini_crop)
        exc_norm_tgw_ini_crop = tgw_SMBH_ini_crop[index_excess_crop] / max(tgw_SMBH_ini_crop)
        mer_norm_tgw_fin_crop = tgw_SMBH_fin_crop[index_merger_crop] / max(tgw_SMBH_ini_crop)
        mer_norm_tgw_ini_crop = tgw_SMBH_ini_crop[index_merger_crop] / max(tgw_SMBH_ini_crop)

        xlist = np.linspace(min(np.log10(tgw_SMBH_ini_crop / max(tgw_SMBH_ini_crop))), max(np.log10(tgw_SMBH_ini / max(tgw_SMBH_ini))))

        cbar_min = min(tcomp_crop)
        cbar_max = max(tcomp_crop)
        normalise = plt.Normalize(np.log10(cbar_min), np.log10(cbar_max))

        plt.figure(figsize=(12, 4))
        gs = gridspec.GridSpec(1, 2)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        plt.set_title('Evolution of Merging Timescale \n for Individual Particles')
        
        ax1.set_ylabel(r'$log_{10}(t_{GW,f}/$max $t_{GW,0})$')
        ax1.set_xlabel(r'$log_{10}(t_{GW,0}/$max $t_{GW,0})$')
        ax2.set_ylabel(r'$log_{10}(t_{GW,f}/t_H)$')
        ax2.set_xlabel(r'$log_{10}(t_{GW,0}/t_H)$')
        for ax_ in [ax1, ax2]:
            plot_init.tickers(ax_)
        ax1.scatter(np.log10(exc_norm_tgw_ini_crop), np.log10(exc_norm_tgw_fin_crop), edgecolors = 'black', c = np.log10(tcomp_crop[index_excess_crop]), norm = normalise)
        ax1.scatter(np.log10(mer_norm_tgw_ini_crop), (np.log10(mer_norm_tgw_fin_crop)), marker = 'X', c = np.log10(tcomp_crop[index_merger_crop]), norm = normalise)
        ax1.plot(xlist, xlist, linestyle = ':', color = 'black')
        plt.colorbar(ax=ax1).set_label(r'$\log_{10}(t_{GW}/t_H)$') #Change to tGW/tbin for later

        xlist = np.linspace(0, 1.1*max(np.log10(tgw_SMBH_ini_ratio)))

        bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(np.log10(tgw_SMBH_fin_ratio), np.log10(tgw_SMBH_ini_ratio), bins=(50,50), \
                                                         range=([0, 1.1*np.log10(max(tgw_SMBH_fin_ratio))],[0, 1.1*np.log10(max(tgw_SMBH_ini_ratio))]))
        bin2d_sim /= np.max(bin2d_sim)
        extent = [0, 1.1*np.log10(max(tgw_SMBH_fin_ratio)), 0, 1.1*np.log10(max(tgw_SMBH_ini_ratio))]
        contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')

        ax2.plot(xlist, xlist, linestyle = ':', color = 'white')
        plt.savefig('figures/tgw_SMBH_evol.pdf', dpi=300, bbox_inches='tight')

        #DO THE SAME FOR GRX
        fig, ax = plt.subplots()
        bin2d_sim, xedges_s, yedges_s, image = ax.hist2d(semi_SMBH, ecc_SMBH, bins=(50,50), range=([0, 1.1*max(semi_SMBH)],[min(ecc_SMBH), 0]))
        bin2d_sim /= np.max(bin2d_sim)
        extent = [-0.1, 1.2*max(semi_SMBH), 1.2*min(ecc_SMBH), -0.1]
        contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
        ax.set_ylabel(r'$\log_{10}(1-e)$')
        ax.set_xlabel(r'$a$ [pc]')
        ax.set_ylim(min(ecc_SMBH), 0)
        plot_init.tickers(ax)
        plt.savefig('figures/ecc_semi_SMBH_histogram.pdf', dpi=300, bbox_inches='tight')

        #DO THE SAME FOR GRX
        fig, ax = plt.subplots()
        bin2d_sim, xedges_s, yedges_s, image = ax.hist2d(np.log10(self.close_enc), np.log10(tgw_evol_SMBH), bins=(50,50), range=([0, 1.1*max(np.log10(self.close_enc))],[min(np.log10(tgw_evol_SMBH)), 1.1*np.max(np.log10(tgw_evol_SMBH))]))
        bin2d_sim /= np.max(bin2d_sim)
        extent = [-0.1, 1.2*max(np.log10(self.close_enc)), 1.2*min(np.log10(tgw_evol_SMBH)), 1.2*np.max(np.log10(tgw_evol_SMBH))]
        contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
        ax.set_xlabel(r'$N_{enc}$')
        ax.set_ylabel(r'$\log_{10}(t_{GW,f}/t_{GW,0})$')
        ax.set_xlim(0, 1.1*max(np.log10(self.close_enc)))
        plot_init.tickers(ax)
        plt.savefig('figures/semi_tgw_Smbh_evolution_histogram.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

        ini_pop = np.unique(self.init_pop)
        avg_tgw = np.empty(len(ini_pop))
        close_enc = np.empty(len(ini_pop))
        iter = -1
        for pop_ in ini_pop:
            iter +=1
            idx = np.where(self.init_pop == pop_)[0]
            temp = [np.log10(tgw_evol_SMBH[i]) for i in idx]
            avg_tgw[iter] = np.mean(temp)
            close_enc[iter] = np.mean(np.asarray(self.close_enc)[idx])

        fig, ax = plt.subplots()
        colour_axes = ax.scatter(ini_pop, avg_tgw, edgecolors = 'black', c = close_enc)
        plot_init.tickers(ax)
        ax.set_xlabel(r'IMBH Population [$N$]')
        ax.set_ylabel(r'$\log_{10}(t_{GW,f}/t_{GW,0}$)')
        ax.set_ylim(0, 1.1*max(avg_tgw))
        plot_init.tickers(ax)
        plt.colorbar(colour_axes, ax = ax, label = r'# Close Encounters')
        plt.savefig('figures/pop_tgwevol_enc_plot.pdf', dpi=300, bbox_inches='tight')

    def semi_ecc_bin_gw_plotter(self):
        """
        Function which detects all binaries and plots their eccentricity vs. semi-major axis
        """

        plot_init = plotter_setup()

        tgw_bin_fin = []
        tgw_bin_SMBH = []

        for i in range(len(self.semi_bin_fin)):
            grav_ftimescale = self.gw_calc(self.semi_bin_fin[i], self.ecc_bin_fin[i], self.mass_bin[i][0], self.mass_bin[i][1]).value_in(units.s)
            grav_SMBH_timescale = self.gw_calc(self.semi_bin_SMBH[i], self.ecc_bin_SMBH[i], self.mass_bin_SMBH[i][0], self.mass_bin_SMBH[i][1]).value_in(units.s)
            tgw_bin_fin.append(grav_ftimescale)
            tgw_bin_SMBH.append(grav_SMBH_timescale)
        tgw_evol_bin = np.asarray([i/j for i, j in zip(tgw_bin_fin, tgw_bin_SMBH)])

        tcomp = np.asarray([i/self.tH for i in min(tgw_bin_fin, tgw_bin_SMBH)])
        semi_bin = np.asarray([i.value_in(units.parsec) for i in self.semi_bin_fin])
        ecc_bin  = np.asarray([np.log10(i) for i in np.asarray(self.ecc_bin_fin)])

        index_excess = np.where(tcomp >= 1)[0]
        index_merger = np.where(tcomp < 1)[0]
        min_norm = min(np.min(tgw_evol_bin[index_excess]), np.min(tgw_evol_bin[index_merger]))
        max_norm = max(np.max(tgw_evol_bin[index_excess]), np.max(tgw_evol_bin[index_merger]))
        normalise = plt.Normalize(np.log10(min_norm), np.log10(max_norm))

        #DO THE SAME FOR GRX
        fig, ax = plt.subplots()
        ax.set_ylabel(r'$\log_{10}(1-e)$')
        ax.set_xlabel(r'$\log_{10}a$ [pc]')
        plot_init.tickers(ax)
        plt.scatter(np.log10(semi_bin[index_excess]), (ecc_bin[index_excess]), edgecolors = 'black', c = (np.log10(tgw_evol_bin[index_excess])), norm = normalise)
        plt.scatter(np.log10(semi_bin[index_merger]), (ecc_bin[index_merger]), marker = 'X', edgecolors='black', c = (np.log10(tgw_evol_bin[index_merger])), norm = normalise)
        plt.colorbar().set_label(r'$\log_{10}(t_{GW,bin}/t_{GW,SMBH})$') #Change to tGW/tbin for later
        plt.savefig('figures/semi_ecc_bin_plot.pdf', dpi=300, bbox_inches='tight')

"""
cst = bin_tert_systems()
cst.semi_ecc_SMBH_gw_plotter()"""