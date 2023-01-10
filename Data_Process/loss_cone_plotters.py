from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pickle as pkl

np.seterr(divide='ignore')
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

class loss_cone(object):
    """
    Class to plot the loss cone evolution
    """

    def __init__(self):
        
        Hermite_data = glob.glob(('/media/erwanh/Elements/Hermite/particle_trajectory/*'))
        GRX_data = glob.glob(('/media/erwanh/Elements/GRX/particle_trajectory/*'))
        pfile = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)]

        cfile = [[], []]
        cfile[0] = ['data/Hermite/no_addition/chaotic_simulation/'+str(i[51:]) for i in pfile[0]]
        cfile[1] = ['data/GRX/no_addition/chaotic_simulation/'+str(i[47:]) for i in pfile[1]]

        dir = os.path.join('/home/erwanh/Desktop/SteadyStateBH/Data_Process/figures/steady_time/Sim_summary.txt')
        with open(dir) as f:
            line = f.readlines()

            popG = line[11][48:-2] 
            avgG = line[17][48:-2]
            avgG2 = line[18][3:-2] 
            popG_data = popG.split()
            avgG_data = avgG.split()
            avgG2_data = avgG2.split()
            popG = np.asarray([float(i) for i in popG_data])
            avgG = np.asarray([float(i) for i in avgG_data])
            avgG2 = np.asarray([float(i) for i in avgG2_data])
            avgG = np.concatenate((avgG, avgG2))

        self.no_data = [len(Hermite_data), len(GRX_data)]
        with open(pfile[0][0], 'rb') as input_file:
            data = pkl.load(input_file)
        self.SMBH_mass = data.iloc[0][1][1]
        self.SMBH_angL = np.sqrt(16*constants.G**2*self.SMBH_mass**2*constants.c**-2).value_in(units.m*units.m/units.s) #TB2008 chapter 7

        self.init_pop = [[ ], [ ]]
        self.angL_init = [[ ], [ ]]
        self.dangL_avg = [[ ], [ ]]
        self.angL_fin = [[ ], [ ]]
        self.KE = [[ ], [ ]]
        self.PE = [[ ], [ ]]
        self.energy = [[ ], [ ]]
        self.angL_evol = [[ ], [ ]]
        self.time = [[ ], [ ]]
        self.ejec_part = [[ ], [ ]]

        for int_ in range(2):
            for file_ in range(len(pfile[int_])):
                with open(cfile[int_][file_], 'rb') as input_file:
                    chaotic_tracker = pkl.load(input_file)
                    pop = chaotic_tracker.iloc[0][6]
                    if pop <= 40 and pop > 5:
                        with open(pfile[int_][file_], 'rb') as input_file:
                            file_size = os.path.getsize(pfile[int_][file_])
                            if file_size < 2.7e9:
                                print('Reading file', file_, ': ', pfile[int_][file_])
                                ptracker = pkl.load(input_file)

                                ###### KEEP THIS FOR COMPARISON PLOTS - FOR EVOLUTION TURN OFF
                                if 10*round(0.1*chaotic_tracker.iloc[0][6]) in popG:
                                    idx = np.where(popG == 10*round(0.1*chaotic_tracker.iloc[0][6]))
                                    col_len = int(min(np.round((avgG[idx])*10**3), np.shape(ptracker)[1])-1)
                                else:
                                    col_len = np.shape(ptracker)[1]-1

                                for parti_ in range(np.shape(ptracker)[0]):
                                    KE_arr = [ ]
                                    PE_arr = [ ]
                                    angL_arr = [ ]
                                    energy_arr = [ ]
                                    time = [ ]
                                    if parti_ != 0:
                                        self.init_pop[int_].append(pop)
                                        angL_val = 0
                                        for col_ in range(col_len):
                                            KE = ptracker.iloc[parti_][col_][4].value_in(units.J)
                                            PE = ptracker.iloc[parti_][col_][5].value_in(units.J)
                                            tot_E = KE + PE
                                            KE_arr.append(KE)
                                            PE_arr.append(PE)
                                            energy_arr.append(abs(tot_E))
                                            time.append(col_*1000)

                                            semi_val = ptracker.iloc[parti_][col_][7][0]
                                            ecc_val = (1-ptracker.iloc[parti_][col_][8][0])
                                            angL_arr.append(self.ang_momentum(semi_val, ecc_val))
                                            angL_val += self.ang_momentum(semi_val, ecc_val)
                                        
                                        if not isinstance(ptracker.iloc[parti_][col_+1][0], np.uint64) or ptracker.iloc[parti_][col_+1][1] > 4*10**6 | units.MSun:
                                            self.ejec_part[int_].append(1)
                                        else:
                                            self.ejec_part[int_].append(0)

                                        sem_init = ptracker.iloc[parti_][0][7][0]
                                        ecc_init = 1-ptracker.iloc[parti_][0][8][0]
                                        semi_fin = ptracker.iloc[parti_][-2][7][0]
                                        ecc_fin = 1-ptracker.iloc[parti_][-2][8][0]

                                        angL_val -= self.ang_momentum(sem_init, ecc_init)
                                        angL_val /= col_len                # Average change per kyr
                                        self.dangL_avg[int_].append(angL_val)

                                        angL_init_val = self.ang_momentum(sem_init, ecc_init)
                                        angL_fin_val = self.ang_momentum(semi_fin, ecc_fin)
                                        self.angL_init[int_].append(angL_init_val)
                                        self.angL_fin[int_].append(angL_fin_val)

                                        self.KE[int_].append(KE_arr)
                                        self.PE[int_].append(PE_arr)
                                        self.energy[int_].append(energy_arr)
                                        self.angL_evol[int_].append(angL_arr)
                                        self.time[int_].append(time)

            self.dangL_avg[int_] = np.asarray(self.dangL_avg[int_])
            self.angL_init[int_] = np.asarray(self.angL_init[int_])
            self.angL_fin[int_] = np.asarray(self.angL_fin[int_]) 
            self.angL_evol[int_] = np.asarray(self.angL_evol[int_])
            self.energy[int_] = np.asarray(self.energy[int_])
            self.ejec_part[int_] = np.asarray(self.ejec_part[int_])
            self.time[int_] = np.asarray(self.time[int_])

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

    def lcone_fininit_plotter(self):
        """
        Plotting function aimed to show the loss-cone evolution during all simulations
        """

        plot_ini = plotter_setup()

        angL_init = np.asarray(self.angL_init)
        angL_fin = np.asarray(self.angL_fin)    
        init_pop_t = np.asarray(self.init_pop)

        init_pop = [[ ], [ ]]
        for int_ in range(2):
            for pop_ in init_pop_t[int_]:
                init_pop[int_].append(10*round(0.1*pop_))
            init_pop[int_] = np.asarray([i for i in init_pop[int_]])
        xline = np.linspace(1, 1.1*max(max(angL_init[0]), max(angL_init[1]),
                                       max(angL_fin[0]),  max(angL_fin[1])))
        
        colours = ['black', 'white']
        file = ['Hermite', 'GRX']
        fig = plt.figure(figsize=(5, 15))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        for int_ in range(2):
            ax_ = [ax1, ax2]
            iter = -1

            xmax = 1.01*np.log10(max(max(angL_init[0]), max(angL_init[1])))
            ymax = 1.01*np.log10(max(max(angL_fin[0]), max(angL_fin[1])))
            plot_ini.tickers(ax1, 'plot')

            bin2d_sim, xed, yed, image = ax_[int_].hist2d(np.log10(angL_init[int_]), 
                                                          np.log10(angL_fin[int_]), bins=(50, 50),
                                                          range=([17.8, 1.005*xmax], [17.8, 1.005*ymax]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [17.8, 1.005*(xmax), 17.8, 1.005*(ymax)]
            contours = ax2.imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')

        for ax_ in [ax1, ax2]:
            plot_ini.tickers(ax_, 'hist')
            ax_.set_xlabel(r'$\log_{10} L_{0}$')
            ax_.set_ylabel(r'$\log_{10} L_{f}$')
            ax_.axhline(np.log10(self.SMBH_angL), color = 'white')
            ax_.axvline(np.log10(self.SMBH_angL), color = 'white')
            ax_.plot(xline, xline, color = 'white', linestyle = ':')
            ax_.set_xlim(17.8, 1.01*(xmax))
            ax_.set_ylim(17.8, 1.01*(ymax))
            ax_.text(18, 1.002*np.log10(self.SMBH_angL), r'$L_{\rm{crit}}$')
        ax1.text(18, 1.002*np.log10(self.SMBH_angL), r'$L_{\rm{crit}}$')

        plt.colorbar(colour_axes, ax=ax1, label = r'IMBH Population [$N$]')
        plt.savefig('figures/loss_cone/LCif_hist_'+str(file[int_])+'.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

        fig = plt.figure(figsize=(5, 15))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        normalise = plt.Normalize(10, 40)
        for int_ in range(2):
            idx_ejec = np.where(self.ejec_part[int_] == 1)[0]
            idx_stab = np.where(self.ejec_part[int_] != 1)[0]
            for idx in idx_ejec:
                colour_axes = ax1.scatter(np.log10(angL_init[int_][idx]), 
                                          np.log10(angL_fin[int_][idx]), 
                                          c = init_pop[int_][idx], norm = normalise,
                                          marker = 'X', edgecolors='black')
            for idx in idx_stab:
                colour_axes = ax1.scatter(np.log10(angL_init[int_][idx]), 
                                          np.log10(angL_fin[int_][idx]), 
                                          c = init_pop[int_][idx], norm = normalise,
                                          edgecolors='black')
            for ax_ in [ax1, ax2]:
                plot_ini.tickers(ax_, 'plot')
                ax_.set_xlabel(r'$\log_{10} L_{0}$')
                ax_.set_ylabel(r'$\log_{10} L_{f}$')
                ax_.axhline(np.log10(self.SMBH_angL), color = colours[iter])
                ax_.axvline(np.log10(self.SMBH_angL), color = colours[iter])
                ax_.plot(xline, xline, color = colours[iter], linestyle = ':')
                ax_.set_xlim(17.8, 1.005*(xmax))
                ax_.set_ylim(17.8, 1.005*(ymax))
            ax1.text(18, 1.002*np.log10(self.SMBH_angL), r'$L_{\rm{crit}}$')

            plt.colorbar(colour_axes, ax=ax1, label = r'IMBH Population [$N$]')
            plt.savefig('figures/loss_cone/LCif_scatter_'+str(file[int_])+'.pdf', dpi=300, bbox_inches='tight')
            plt.clf()

    def lcone_plotter(self):
        """
        Function to plot the change in energy and angular momentum of agiven particle
        """

        plot_ini = plotter_setup()

        fig = plt.figure(figsize=(16, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax = [ax1, ax2]
        normalise = plt.Normalize(0, 8)

        for int_ in range(2):
            idx_ejec = np.where(self.ejec_part[int_] == 1)[0]
            for idx in idx_ejec:
                ax[int_].scatter(np.log10(self.angL_evol[int_][idx][0]), np.log10((self.energy[int_][idx][0])), color = 'black', s = 18, zorder = 4)
                ax[int_].scatter(np.log10(self.angL_evol[int_][idx][-1]), np.log10((self.energy[int_][idx][-1])), color = 'black', s = 18, zorder = 3)
                colour_axes = ax[int_].scatter(np.log10(self.angL_evol[int_][idx][1:-3]), np.log10((self.energy[int_][idx][1:-3])), 
                                                c = np.log10(self.time[int_][idx][1:-3]), norm = normalise, s = 13, zorder = 1)
        for int_ in range(2):
            idx_ejec = np.where(self.ejec_part[int_] == 1)[0]
            for idx in idx_ejec:
                ax[int_].plot(np.log10(self.angL_evol[int_][idx]), np.log10((self.energy[int_][idx])), color = 'black', alpha = 0.4, zorder = 2)
                
        for ax_ in [ax1, ax2]:
            ax_.set_xlabel(r'$L$ [m$^{2}$ s$^{-1}$]')
            ax_.set_ylabel(r'$\log_{10}|E_{\rm{tot}}|$ [J]')
            ax_.set_ylim(43.25, 44.5)
            ax_.set_xlim(17.8, 23)
            plot_ini.tickers(ax_, 'plot')
            ax_.axhline(np.log10(self.SMBH_angL), color = 'black')
            ax_.axvline(np.log10(self.SMBH_angL), color = 'black')
        ax1.text(18.85, 43.4, r'$L_{\rm{crit}}$', rotation = -90)
        
        plt.colorbar(colour_axes, ax = ax2, label = r'$\log_{10} t$ [yr]')
        plt.savefig('figures/loss_cone/LC_evol.png', dpi=300, bbox_inches='tight')

    def lcone_timescale(self):
        init_pop_t = np.asarray(self.init_pop)
        init_pop = [[ ], [ ]]
        for int_ in range(2):
            for pop_ in init_pop_t[int_]:
                init_pop[int_].append(10*round(0.1*pop_))
            init_pop[int_] = np.asarray([i for i in init_pop[int_]])
            init_pop[int_] = np.unique(init_pop[int_])
        avg_angL_timescale = [[ ], [ ]]
        avg_angL_replenish = [[ ], [ ]]
        
        iter = -1
        for int_ in range(2):
            for pop_ in init_pop[int_]:
                iter += 1            
                indices = np.where((init_pop[int_] == pop_))[0]
                tscale = (1000)*np.mean((self.angL_init[int_][indices]-self.SMBH_angL)/self.dangL_avg[int_][indices])
                avg_angL_timescale[int_].append(tscale)
                avg_angL_replenish[int_].append(tscale/pop_)

        plot_ini = plotter_setup()
        colours = ['red', 'blue']
        labels = ['Hermite', 'GRX']
        fig = plt.figure(figsize=(12.5, 5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        plot_ini.tickers_pop(ax1, init_pop[0], 'Hermite')
        plot_ini.tickers_pop(ax2, init_pop[1], 'GRX')
        ax1.set_ylabel(r'$\tau_{\rm{repl}}$ [Myr]')
        ax2.set_ylabel(r'$\langle t_{\rm{LC, repl}}\rangle$ [Myr]')
        for int_ in range(2):
            ax1.scatter(init_pop[int_], avg_angL_replenish[int_], \
                        color = colours[int_], edgecolors = 'black', label = labels[int_])
            ax2.scatter(init_pop[int_], avg_angL_timescale[int_], \
                        color = colours[int_], edgecolors = 'black', label = labels[int_])
            ax1.legend()
        plt.savefig('figures/loss_cone/avgN_LC_tscale_repl.pdf', dpi=300, bbox_inches='tight')