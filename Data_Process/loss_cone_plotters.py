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
        self.IMBH_tracker = bulk_stat_extractor('data/Hermite/particle_trajectory_temp/*', 'Y')
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

        for j in range(self.no_data):
            sim_data = self.IMBH_tracker[j]
            for parti_ in range(len(sim_data)):
                if parti_ == 0:  # Don't care for SMBH
                    pass
                else:
                    self.init_pop.append((np.shape(self.IMBH_tracker[j]))[0])
                    if isinstance(sim_data.iloc[parti_][-1][0], np.uint64):           # Neglect removed particle
                        angL_val = 0
                        for col_ in range(np.shape(sim_data)[1]-1):
                            semi_val = sim_data.iloc[parti_][col_][7][0]
                            ecc_val = (1-sim_data.iloc[parti_][col_][8][0])
                            angL_val += self.ang_momentum(semi_val, ecc_val)
                        angL_val /= np.shape(sim_data)[1]
                        self.angL_avg.append(angL_val)

                        semi_val_init = sim_data.iloc[parti_][2][7][0]
                        ecc_val_init = 1-sim_data.iloc[parti_][2][8][0]
                        angL_init_val = self.ang_momentum(semi_val_init, ecc_val_init)
                        self.angL_init.append(angL_init_val)

                        semi_val_fin = sim_data.iloc[parti_][-2][7][0]
                        ecc_val_fin = 1-sim_data.iloc[parti_][-2][8][0]
                        angL_fin_val = self.ang_momentum(semi_val_fin, ecc_val_fin)
                        self.angL_fin.append(angL_fin_val)
                    else:
                        angL_val = 0
                        for col_ in range(np.shape(sim_data)[1]-1):
                            semi_val = sim_data.iloc[parti_][col_][7][0]
                            ecc_val = (1-sim_data.iloc[parti_][col_][8][0])
                            angL_val += self.ang_momentum(semi_val, ecc_val)
                        angL_val /= np.shape(sim_data)[1]
                        self.angL_avg.append(angL_val)

                        semi_val_init = sim_data.iloc[parti_][2][7][0]
                        ecc_val_init = 1-sim_data.iloc[parti_][2][8][0]
                        angL_init_val = self.ang_momentum(semi_val_init, ecc_val_init)
                        self.angL_init.append(angL_init_val)

                        semi_val_fin = sim_data.iloc[parti_][-3][7][0]
                        ecc_val_fin = 1-sim_data.iloc[parti_][-3][8][0]
                        angL_fin_val = self.ang_momentum(semi_val_fin, ecc_val_fin)
                        self.angL_fin.append(angL_fin_val)

        return self.angL_avg, self.angL_init, self.angL_fin

    def log_ratio(self):
        angL_avg, angL_init, angL_fin = self.data_extractor()

        angL_avg = np.log10(self.angL_avg)
        angL_init = np.log10(self.angL_init)
        angL_init /= self.SMBH_angL
        angL_fin = np.log10(self.angL_fin)
        angL_fin /= self.SMBH_angL

        return angL_avg, angL_init, angL_fin

    def lcone_fininit_plotter(self):

        plot_ini = plotter_setup()

        angL_avg, angL_init, angL_fin = self.log_ratio()
        angL_avg = np.asarray(angL_avg)
        
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
            plot_ini.tickers(ax_)

        bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(angL_init, angL_fin, bins=(20, 20), range=([0.75*min(angL_init),1.1*max(angL_init)],[0.75*min(angL_fin), 1.1*max(angL_fin)]))
        extent = [0,max(angL_init), 0, max(angL_fin)]
        bin2d_sim /= np.max(bin2d_sim)
        contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
        colour_axes = ax1.scatter(angL_init, angL_fin, c = self.init_pop, edgecolors='black')
        plt.colorbar(colour_axes, ax=ax1, label = r'IMBH Population [$N$]')
        plt.savefig('figures/loss_cone_evolution.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

    def lcone_evolution_plotter(self):

        angL_avg, angL_init, angL_fin = self.log_ratio
        angL_avg = np.asarray(angL_avg)

        init_pop = []
        for pop_ in self.init_pop:
            if pop_%10 == 0:
                init_pop.append(pop_)
            else:
                init_pop.append(pop_-1)
        init_pop = np.unique(init_pop)
        avg_angL = np.empty(len(init_pop))

        iter = -1
        for pop_ in init_pop:
            iter += 1            
            indices = np.where((init_pop == pop_))[0]
            avg_angL[iter] = np.mean(angL_avg[indices])
        avg_angL /= max(avg_angL)

        plot_ini = plotter_setup()
        fig, ax = plt.subplots()
        plot_ini.tickers_pop(ax, self.init_pop)
        ax.set_ylabel(r'$\langle \Delta L\rangle / \langle \Delta L\rangle_{max}$')
        ax.scatter(init_pop, avg_angL, color = 'red', edgecolors = 'black', label = 'Hermite')
        ax.legend()
        plt.savefig('figures/avg_N_loss_cone_evolution.pdf', dpi=300, bbox_inches='tight')

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

        iter = -1
        for pop_ in init_pop:
            iter += 1            
            indices = np.where((init_pop == pop_))[0]
            avg_angL_timescale[iter] = 1000*np.mean(angL_init[indices]-self.SMBH_angL)/np.mean(angL_avg[indices])

        plot_ini = plotter_setup()
        fig, ax = plt.subplots()
        plot_ini.tickers_pop(ax, self.init_pop)
        ax.set_ylabel(r'$t_{LC, repl}$ [Myr]')
        ax.scatter(init_pop, avg_angL_timescale, color = 'red', edgecolors = 'black', label = 'Hermite')
        #ax.hline()
        ax.legend()
        plt.savefig('figures/avg_N_loss_cone_evolution.pdf', dpi=300, bbox_inches='tight')

        
cst = loss_cone()
cst.lcone_timescale()
cst.lcone_fininit_plotter()
cst.lcone_evolution_plotter()