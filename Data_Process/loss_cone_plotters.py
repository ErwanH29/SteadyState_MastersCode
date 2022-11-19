from amuse.lab import *
from spatial_plotters import *
from file_logistics import *
import matplotlib.pyplot as plt
import numpy as np
import warnings

np.seterr(divide='ignore')
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

class loss_cone(object):
    """
    Class to plot the loss cone evolution
    """

    def __init__(self):
        self.IMBH_tracker = bulk_stat_extractor('data/Hermite/particle_trajectory/*', 'Y')
        self.IMBH_tracker_GRX = bulk_stat_extractor('data/GRX/particle_trajectory/*', 'Y')
        self.no_data = [len(self.IMBH_tracker), len(self.IMBH_tracker_GRX)]
        self.data_array = [self.IMBH_tracker, self.IMBH_tracker_GRX]

        self.SMBH_mass = self.IMBH_tracker[0].iloc[0][1][1]
        self.SMBH_angL = np.log10(np.sqrt((((2*constants.G*self.SMBH_mass)*(8*constants.G*self.SMBH_mass*constants.c**-2)))).value_in(units.m*units.m/units.s)) #TB2008 chapter 7

    def ang_momentum(self, semi_ax, ecc):
        """
        Function to calculate the angular momentum
        
        Inputs:
        semi_ax: The particles semi-axis with the SMBH at any given time
        ecc:     The particles eccentricity with the SMBH at any given time
        outputs: The angular momentum of a particular particle
        """
        qK_val = semi_ax * ecc
        angMom_val = (2*constants.G*self.SMBH_mass*qK_val).sqrt().value_in(units.m*units.m/units.s) #TB2008 chapter 7
        return angMom_val

    def data_extractor(self):
        """
        Function which extracts needed data for the angular momentum calculations
        """
        self.init_pop = [[ ], [ ]]
        self.angL_init = [[ ], [ ]]
        self.angL_avg = [[ ], [ ]]
        self.angL_fin = [[ ], [ ]]

        for int_ in range(2):
            for j in range(self.no_data[int_]):
                sim_data = self.data_array[int_][j]
                for parti_ in range(len(sim_data)):
                    if parti_ !=0 and isinstance(sim_data.iloc[parti_][-1][0], np.uint64):
                        pop = 10*round(0.1*(np.shape(self.data_array[int_][j]))[0])
                        self.init_pop[int_].append(pop)
                        angL_val = 0
                        for col_ in range(np.shape(sim_data)[1]-1):
                            semi_val = sim_data.iloc[parti_][col_][7][0]
                            ecc_val = (1-sim_data.iloc[parti_][col_][8][0])
                            angL_val += self.ang_momentum(semi_val, ecc_val)
                        angL_val /= np.shape(sim_data)[1] - 1                # Average change per kyr
                        self.angL_avg[int_].append(angL_val)

                        semi_val_init = sim_data.iloc[parti_][2][7][0]
                        ecc_val_init = 1-sim_data.iloc[parti_][2][8][0]
                        angL_init_val = self.ang_momentum(semi_val_init, ecc_val_init)
                        self.angL_init[int_].append(angL_init_val)

                        semi_val_fin = sim_data.iloc[parti_][-2][7][0]
                        ecc_val_fin = 1-sim_data.iloc[parti_][-2][8][0]
                        angL_fin_val = self.ang_momentum(semi_val_fin, ecc_val_fin)
                        self.angL_fin[int_].append(angL_fin_val)
            self.angL_avg[int_] = np.asarray(self.angL_avg[int_]) 
            self.angL_init[int_] = np.asarray(self.angL_init[int_]) 
            self.angL_fin[int_] = np.asarray(self.angL_fin[int_]) 
            self.init_pop[int_] = np.asarray(self.init_pop[int_])
                            
        return self.angL_avg, self.angL_init, self.angL_fin

    def log_ratio(self):
        """
        Function to normalise all the angular momentum data to make it more presentable
        """

        angL_avg, angL_init, angL_fin = self.data_extractor()

        for int_ in range(2):
            angL_avg[int_] = np.log10(self.angL_avg[int_])
            angL_init[int_] = np.log10(self.angL_init[int_])
            angL_fin[int_] = np.log10(self.angL_fin[int_])
            angL_init[int_] /= self.SMBH_angL
            angL_fin[int_] /= self.SMBH_angL

        return angL_avg, angL_init, angL_fin

    def lcone_fininit_plotter(self):
        """
        Plotting function aimed to show the loss-cone evolution during all simulations
        """

        plot_ini = plotter_setup()

        angL_avg, angL_init, angL_fin = self.log_ratio()
        angL_avg = np.asarray(angL_avg)
        angL_init = np.asarray(angL_init)
        angL_fin = np.asarray(angL_fin)
        init_pop = np.asarray(self.init_pop)

        xline = np.linspace(1, 1.1*max(max(angL_init[0]), max(angL_init[1]),
                                       max(angL_fin[0]), max(angL_fin[1])))
        file = ['Hermite', 'GRX']
        for int_ in range(2):
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            ax_ = [ax1, ax2]

            colours = ['black', 'white']
            xmin = 0.95*min(min(angL_init[0]), min(angL_init[1]))
            xmax = 1.01*max(max(angL_init[0]), max(angL_init[1]))
            ymin = 0.95*min(min(angL_fin[0]), min(angL_fin[1]))
            ymax = 1.01*max(max(angL_fin[0]), max(angL_fin[1]))
            plot_ini.tickers(ax2, 'hist')
            plot_ini.tickers(ax1, 'plot')

            angL_itemp = np.asarray([i**10 for i in angL_init[int_]])
            angL_ftemp = np.asarray([i**10 for i in angL_fin[int_]])
            filter_val = (abs(angL_itemp - angL_ftemp)/angL_itemp)
            idx = np.argwhere(filter_val > 5e-1)
        
            bin2d_sim, xed, yed, image = ax2.hist2d(angL_init[int_], angL_fin[int_], bins=(100, 100), \
                                                    range=([0.98*xmin, 1.02*xmax], [0.98*ymin, 1.02*ymax]))
            extent = [0.98*xmin, 1.02*xmax, 0.98*ymin, 1.02*ymax]
            bin2d_sim /= np.max(bin2d_sim)
            contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            colour_axes = ax1.scatter(angL_init[int_][idx], angL_fin[int_][idx], c = init_pop[int_][idx], edgecolors='black')
            ax1.text(1.01, 1.01, r'$\|\Delta L\|/L_i > 0.5$')

            iter = -1
            for ax_ in [ax1, ax2]:
                iter += 1
                ax_.set_xlabel(r'$\log_{10}(L_{0} / L_{\rm{crit}})$')
                ax_.set_ylabel(r'$\log_{10}(L_{f} / L_{\rm{crit}})$')
                ax_.axhline(1, color = colours[iter])
                ax_.axvline(1, color = colours[iter])
                ax_.plot(xline, xline, color = colours[iter], linestyle = ':')
                ax_.set_xlim(0.98*xmin, 1.02*xmax)
                ax_.set_ylim(0.98*ymin, 1.02*ymax)

            plt.colorbar(colour_axes, ax=ax1, label = r'IMBH Population [$N$]')
            plt.savefig('figures/loss_cone/loss_cone_evolution'+str(file[int_])+'.pdf', dpi=300, bbox_inches='tight')
            plt.clf()

    def lcone_timescale(self):
        angL_avg, angL_init, angL_fin = self.data_extractor()
        init_pop_t = np.asarray(self.init_pop)

        init_pop = [[ ], [ ]]
        for int_ in range(2):
            for pop_ in init_pop_t[int_]:
                if pop_%10 == 0:
                    init_pop[int_].append(pop_)
                else:
                    init_pop[int_].append(pop_-1)
            init_pop[int_] = np.unique(init_pop[int_])
        avg_angL_timescale = [[ ], [ ]]
        avg_angL_replenish = [[ ], [ ]]
        
        iter = -1
        for int_ in range(2):
            for pop_ in init_pop[int_]:
                iter += 1            
                indices = np.where((init_pop[int_] == pop_))[0]
                tscale = (1000)*np.mean((angL_init[int_][indices]-self.SMBH_angL)/angL_avg[int_][indices])
                avg_angL_timescale[int_].append(tscale)
                avg_angL_replenish[int_].append(tscale/pop_)

        plot_ini = plotter_setup()
        colours = ['red', 'blue']
        labels = ['Hermite', 'GRX']
        fig = plt.figure(figsize=(12.5, 5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        for ax_ in [ax1, ax2]:
            plot_ini.tickers_pop(ax_, init_pop[0])
        ax1.set_ylabel(r'$\tau_{\rm{repl}}$ [Myr]')
        ax2.set_ylabel(r'$\langle t_{\rm{LC, repl}}\rangle$ [Myr]')
        for int_ in range(2):
            ax1.scatter(init_pop[int_], avg_angL_replenish[int_], \
                        color = colours[int_], edgecolors = 'black', label = labels[int_])
            ax2.scatter(init_pop[int_], avg_angL_timescale[int_], \
                        color = colours[int_], edgecolors = 'black', label = labels[int_])
            ax1.legend()
        plt.savefig('figures/loss_cone/avg_N_loss_cone_tscale_repl.pdf', dpi=300, bbox_inches='tight')

"""print('...loss_cone_plotters...')
cst = loss_cone()
cst.lcone_timescale()
cst.lcone_fininit_plotter()"""