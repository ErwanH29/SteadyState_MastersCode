from amuse.lab import *
from spatial_plotters import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

class loss_cone(object):
    """
    Class to plot the loss cone evolution
    """

    def __init__(self):
        self.IMBH_tracker = bulk_stat_extractor('data/Hermite/particle_trajectory/*', 'Y')
        self.no_data = len(self.IMBH_tracker)

        self.SMBH_mass = self.IMBH_tracker[0].iloc[0][1][1]
        self.SMBH_angL = np.log10(np.sqrt((((2*constants.G*self.SMBH_mass)*(8*constants.G*self.SMBH_mass*constants.c**-2)))).value_in(units.m*units.m/units.s))

        self.angL_init = []
        self.angL_avg = []
        self.angL_fin = []
        self.init_pop = []

    def ang_momentum(self, semi_ax, ecc):
        """
        Function to calculate the angular momentum
        
        Inputs:
        semi_ax: The particles semi-axis with the SMBH at any given time
        ecc:     The particles eccentricity with the SMBH at any given time
        outputs: The angular momentum of a particular particle
        """
        qK_val = semi_ax * ecc
        angMom_val = (2*constants.G*self.SMBH_mass*qK_val).sqrt().value_in(units.m*units.m/units.s)
        return angMom_val

    def data_extractor(self):
        """
        Function which extracts needed data for the angular momentum calculations
        """

        for j in range(self.no_data):
            sim_data = self.IMBH_tracker[j]
            for parti_ in range(len(sim_data)):
                if parti_ == 0:  # Don't care for SMBH
                    pass
                else:
                    if isinstance(sim_data.iloc[parti_][-1][0], np.uint64):           # Neglect removed particle
                        self.init_pop.append((np.shape(self.IMBH_tracker[j]))[0])
                        angL_val = 0
                        for col_ in range(np.shape(sim_data)[1]-1):
                            semi_val = sim_data.iloc[parti_][col_][7][0]
                            ecc_val = (1-sim_data.iloc[parti_][col_][8][0])
                            angL_val += self.ang_momentum(semi_val, ecc_val)
                        angL_val /= np.shape(sim_data)[1] - 1                         # Average change per kyr
                        self.angL_avg.append(angL_val)

                        semi_val_init = sim_data.iloc[parti_][2][7][0]
                        ecc_val_init = 1-sim_data.iloc[parti_][2][8][0]
                        angL_init_val = self.ang_momentum(semi_val_init, ecc_val_init)
                        self.angL_init.append(angL_init_val)

                        semi_val_fin = sim_data.iloc[parti_][-2][7][0]
                        ecc_val_fin = 1-sim_data.iloc[parti_][-2][8][0]
                        angL_fin_val = self.ang_momentum(semi_val_fin, ecc_val_fin)
                        self.angL_fin.append(angL_fin_val)

        return self.angL_avg, self.angL_init, self.angL_fin

    def log_ratio(self):
        """
        Function to normalise all the angular momentum data to make it more presentable
        """

        angL_avg, angL_init, angL_fin = self.data_extractor()

        angL_avg = np.log10(self.angL_avg)
        angL_init = np.log10(self.angL_init)
        angL_init /= self.SMBH_angL
        angL_fin = np.log10(self.angL_fin)
        angL_fin /= self.SMBH_angL

        return angL_avg, angL_init, angL_fin

    def lcone_fininit_plotter(self):
        """
        Plotting function aimed to show the loss-cone evolution during all simulations
        """

        plot_ini = plotter_setup()

        angL_avg, angL_init, angL_fin = self.log_ratio()
        angL_avg = np.asarray(angL_avg)
        angL_init = np.asarray(angL_init)
        init_pop = np.asarray(self.init_pop)

        angL_itemp = np.asarray([i**10 for i in angL_init])
        angL_ftemp = np.asarray([i**10 for i in angL_fin])
        filter_val = (abs(angL_itemp - angL_ftemp)/angL_itemp)
        idx = np.argwhere(filter_val > 5e-1)
        
        xline = np.linspace(1, 1.1*max(angL_init))

        plt.figure(figsize=(12, 4))
        gs = gridspec.GridSpec(1, 2)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        colours = ['black', 'white']
        iter = -1
        for ax_ in [ax1, ax2]:
            iter += 1
            ax_.set_xlabel(r'$\log_{10}(L_{0} / L_{crit})$')
            ax_.set_ylabel(r'$\log_{10}(L_{f} / L_{crit})$')
            ax_.axhline(1, color = colours[iter])
            ax_.axvline(1, color = colours[iter])
            ax_.plot(xline, xline, color = colours[iter], linestyle = ':')
            ax_.set_xlim(0.95*min(angL_init), 1.01*max(angL_init))
            ax_.set_ylim(0.95*min(angL_fin),  1.01*max(angL_fin))
        plot_ini.tickers(ax2, 'hist')
        plot_ini.tickers(ax2, 'plot')

        bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(angL_init, angL_fin, bins=(100, 100), range=([0.95*min(angL_init),1.01*max(angL_init)],[0.95*min(angL_fin), 1.01*max(angL_fin)]))
        extent = [0,1.01*max(angL_init), 0, 1.01*max(angL_fin)]
        bin2d_sim /= np.max(bin2d_sim)
        contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
        colour_axes = ax1.scatter(angL_init[idx], angL_fin[idx], c = init_pop[idx], edgecolors='black')
        ax1.text(1.1, 1.01, r'$\|\Delta L\|/L_i > 0.5$')
        plt.colorbar(colour_axes, ax=ax1, label = r'IMBH Population [$N$]')
        plt.savefig('figures/loss_cone_evolution.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

    def lcone_timescale(self):
        angL_avg, angL_init, angL_fin = self.data_extractor()
        angL_avg = np.asarray(angL_avg)
        angL_init = np.asarray(angL_init)
        angL_fin = np.asarray(angL_fin)

        init_pop = []
        for pop_ in self.init_pop:
            if pop_%10 == 0:
                init_pop.append(pop_)
            else:
                init_pop.append(pop_-1)
        init_pop = np.unique(init_pop)
        avg_angL_timescale = np.empty(len(init_pop))
        avg_angL_replenish = np.empty(len(init_pop))

        iter = -1
        for pop_ in init_pop:
            iter += 1            
            indices = np.where((init_pop == pop_))[0]
            avg_angL_timescale[iter] = (1000)*np.mean(angL_init[indices]-self.SMBH_angL)/np.mean(angL_avg[indices])
            avg_angL_replenish[iter] = avg_angL_timescale[iter]/pop_

        plot_ini = plotter_setup()
        fig = plt.figure(figsize=(12.5, 5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        for ax_ in [ax1, ax2]:
            plot_ini.tickers_pop(ax_, self.init_pop)
        ax1.set_ylabel(r'$t_{repl,rate}$ [Myr]')
        ax1.scatter(init_pop, avg_angL_replenish, color = 'red', edgecolors = 'black', label = 'Hermite')
        ax2.set_ylabel(r'$\langle t_{LC, repl}\rangle$ [Myr]')
        ax2.scatter(init_pop, avg_angL_timescale, color = 'red', edgecolors = 'black')
        ax1.legend()
        plt.savefig('figures/avg_N_loss_cone_tscale_repl.pdf', dpi=300, bbox_inches='tight')
