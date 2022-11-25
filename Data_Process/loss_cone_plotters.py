from amuse.lab import *
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
        
        Hermite_data = glob.glob('data/Hermite/particle_trajectory/*')
        GRX_data = glob.glob('data/GRX/particle_trajectory/*')
        filename = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)] 
        GRX_chaotic = ['data/GRX/no_addition/chaotic_simulation/'+str(i[29:]) for i in filename[1]]
        chaotic = glob.glob('data/Hermite/no_addition/chaotic_simulation/*')
        filename_c = [natsort.natsorted(chaotic), natsort.natsorted(GRX_chaotic)] 
        self.no_data = [len(Hermite_data), len(GRX_data)]

        with open(filename[0][0], 'rb') as input_file:
            data = pkl.load(input_file)
            self.SMBH_mass = data.iloc[0][1][1]
        self.SMBH_angL = (np.sqrt((((2*constants.G*self.SMBH_mass)*(8*constants.G*self.SMBH_mass*constants.c**-2)))).value_in(units.m*units.m/units.s)) #TB2008 chapter 7

        self.init_pop = [[ ], [ ]]
        self.angL_init = [[ ], [ ]]
        self.angL_avg = [[ ], [ ]]
        self.angL_fin = [[ ], [ ]]
        self.KE = [[ ], [ ]]
        self.PE = [[ ], [ ]]
        self.energy = [[ ], [ ]]
        self.angL_evol = [[ ], [ ]]
        self.time = [[ ], [ ]]
        self.ejected = [[ ], [ ]]

        integrator = ['Hermite', 'GRX']
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename_c[int_][file_], 'rb') as input_file:
                    chaotic_tracker = pkl.load(input_file)
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ': ', filename[int_][file_])
                    data = pkl.load(input_file)
                    pop = 10*round(0.1*(np.shape(data))[0])
                    focus_idx, res = ejected_index(data, chaotic_tracker, integrator[int_])
                    self.ejected[int_].append(focus_idx)
                    for parti_ in range(np.shape(data)[0]):
                        KE_arr = [ ]
                        PE_arr = [ ]
                        angL_arr = [ ]
                        energy_arr = [ ]
                        time = [ ]
                        if parti_ !=0:
                            self.init_pop[int_].append(pop)
                            angL_val = 0
                            for col_ in range(np.shape(data)[1]-1):
                                KE = data.iloc[parti_][col_][4].value_in(units.J)
                                PE = data.iloc[parti_][col_][5].value_in(units.J)
                                tot_E = KE + PE
                                KE_arr.append(KE)
                                PE_arr.append(PE)
                                energy_arr.append(tot_E)
                                time.append(col_*1000)

                                semi_val = data.iloc[parti_][col_][7][0]
                                ecc_val = (1-data.iloc[parti_][col_][8][0])
                                angL_arr.append(self.ang_momentum(semi_val, ecc_val))
                                angL_val += self.ang_momentum(semi_val, ecc_val)
                                
                            angL_val /= np.shape(data)[1] - 1                # Average change per kyr
                            self.angL_avg[int_].append(angL_val)

                            semi_val_init = data.iloc[parti_][2][7][0]
                            ecc_val_init = 1-data.iloc[parti_][2][8][0]
                            angL_init_val = self.ang_momentum(semi_val_init, ecc_val_init)
                            self.angL_init[int_].append(angL_init_val)

                            semi_val_fin = data.iloc[parti_][-2][7][0]
                            ecc_val_fin = 1-data.iloc[parti_][-2][8][0]
                            angL_fin_val = self.ang_momentum(semi_val_fin, ecc_val_fin)
                            self.angL_fin[int_].append(angL_fin_val)

                            self.KE[int_].append(KE_arr)
                            self.PE[int_].append(PE_arr)
                            self.energy[int_].append(energy_arr)
                            self.angL_evol[int_].append(angL_arr)
                            self.time[int_].append(time)

            self.angL_avg[int_] = np.asarray([i for i in self.angL_avg[int_]])
            self.angL_init[int_] = np.asarray([i for i in self.angL_init[int_]])
            self.angL_fin[int_] = np.asarray([i for i in self.angL_fin[int_]]) 
            self.angL_evol[int_] = np.asarray([i for i in self.angL_evol[int_]])
            self.energy[int_] = np.asarray([i for i in self.energy[int_]])

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

    def log_ratio(self):
        """
        Function to normalise all the angular momentum data to make it more presentable
        """

        for int_ in range(2):
            self.angL_avg[int_] = (self.angL_avg[int_])
            self.angL_init[int_] /= self.SMBH_angL
            self.angL_fin[int_] /= self.SMBH_angL
            self.angL_init[int_] = np.asarray([np.log10(i) for i in self.angL_init[int_]])
            self.angL_fin[int_] = np.asarray([np.log10(i) for i in self.angL_fin[int_]])

        return self.angL_avg, self.angL_init, self.angL_fin

    def lcone_fininit_plotter(self):
        """
        Plotting function aimed to show the loss-cone evolution during all simulations
        """

        plot_ini = plotter_setup()

        angL_avg, angL_init, angL_fin = self.log_ratio()
        angL_avg = np.asarray(angL_avg)
        angL_init = np.asarray(angL_init)
        angL_fin = np.asarray(angL_fin)        
        init_pop_t = np.asarray(self.init_pop)
        init_pop = [[ ], [ ]]
        for int_ in range(2):
            for pop_ in init_pop_t[int_]:
                if pop_%10 == 0:
                    init_pop[int_].append(pop_)
                else:
                    init_pop[int_].append(pop_-1)
            init_pop[int_] = np.asarray([i for i in init_pop[int_]])

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
            colour_axes = ax1.scatter(angL_init[int_][idx], 
                                      angL_fin[int_][idx], 
                                      c = init_pop[int_][idx], edgecolors='black')
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

    def lcone_plotter(self):
        """
        Function to plot the change in energy and angular momentum of agiven particle
        """

        plot_ini = plotter_setup()

        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax = [ax1, ax2]
        self.angL_evol /= self.SMBH_angL
        normalise = plt.Normalize(0, 6)

        for ax_ in [ax1, ax2]:
            ax_.set_xlabel(r'$\log_{10}(L_c / L_{LC})$ [m$^{2}$ s$^{-1}$]')
            ax_.set_ylabel(r'$\log_{10}|E_{\rm{tot}}|$ [J]')
            ax_.set_ylim(43.25, 44.5)
            ax_.set_xlim(0,2.8)
            plot_ini.tickers(ax_, 'plot')
        for int_ in range(2):
            for parti_ in self.ejected[int_]:
                ax[int_].plot(np.log10(self.angL_evol[int_][parti_][0:-1]), np.log10(abs(self.energy[int_][parti_][0:-1])), c = 'black', linewidth = 2, alpha = 0.3, zorder = 1)
        for int_ in range(2):
            for parti_ in self.ejected[int_]:
                ax[int_].scatter(np.log10(self.angL_evol[int_][parti_][0]), np.log10(abs(self.energy[int_][parti_][0])), color = 'black', s = 18, zorder = 4)
                ax[int_].scatter(np.log10(self.angL_evol[int_][parti_][-1]), np.log10(abs(self.energy[int_][parti_][-1])), color = 'black', s = 18, zorder = 3)
                colour_axes = ax[int_].scatter(np.log10(self.angL_evol[int_][parti_][1:-3]), np.log10(abs(self.energy[int_][parti_][1:-3])), 
                                               c = np.log10(self.time[int_][parti_][1:-3]), norm = normalise, s = 13, zorder = 2)
                ax[int_].axvline(1, color = 'black')
        plt.colorbar(colour_axes, ax = ax2, label = r'$\log_{10} t$ [Myr]')
        plt.savefig('figures/loss_cone/merger_evol.pdf', dpi=300, bbox_inches='tight')

    def lcone_timescale(self):
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
                tscale = (1000)*np.mean((self.angL_init[int_][indices]-self.SMBH_angL)/self.angL_avg[int_][indices])
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
