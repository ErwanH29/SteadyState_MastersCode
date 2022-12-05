from file_logistics import *
from amuse.ext.galactic_potentials import MWpotentialBovy2015
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import warnings

class ejection_stats(object):
    """
    Class which extracts EJECTED data files
    """
    def __init__(self):
        self.integrator = ['Hermite', 'GRX']
        warnings.filterwarnings("ignore", category=RuntimeWarning) 
        warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

    def new_data_extractor(self):
        """
        Function to extract newly simulated data into reduced manipulatable files
        """

        Hermite_data = glob.glob(os.path.join('/media/erwanh/Elements/Hermite/particle_trajectory_ejec/*'))
        chaotic_H = ['data/Hermite/no_addition/chaotic_simulation/'+str(i[51:]) for i in Hermite_data]
        GRX_data = glob.glob(os.path.join('/media/erwanh/Elements/GRX/particle_trajectory_ejec/*'))
        chaotic_G = ['data/GRX/no_addition/chaotic_simulation/'+str(i[47:]) for i in GRX_data]

        filest = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)] 
        filesc = [natsort.natsorted(chaotic_H), natsort.natsorted(chaotic_G)] 
        for int_ in range(2):
            for file_ in range(len(filest[int_])):
                with open(filesc[int_][file_], 'rb') as input_file:
                    ctracker = pkl.load(input_file)
                    if ctracker.iloc[0][5] > 0:
                        with open(filest[int_][file_], 'rb') as input_file:
                            print('Reading File :', input_file)
                            count = len(fnmatch.filter(os.listdir('data/ejection_stats/'), '*.*'))
                            ptracker = pkl.load(input_file)

                            ex, ey, ez, vesc, ejec_KE, ejec_PE, Nclose, ebool = ejected_extract_final(ptracker, ctracker, 'E')

                            path = '/home/erwanh/Desktop/SteadyStateBH/Data_Process/data/ejection_stats/'
                            stab_tracker = pd.DataFrame()
                            df_stabtime = pd.Series({'self.integrator': self.integrator[int_],
                                                     'Population': np.shape(ptracker)[0],
                                                     'Simulation Time': np.shape(ptracker)[1],
                                                     'xPos': ex, 'yPos': ey, 'zPos': ez,
                                                     'vesc': vesc, 'KE': ejec_KE, 'PE': ejec_PE,
                                                     'Nclose': Nclose, 'Ejected Bool': ebool})
                            stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
                            stab_tracker.to_pickle(os.path.join(path, 'IMBH_'+str(self.integrator[int_])+'_ejec_data_indiv_parti_'+str(count)+'.pkl'))

    def combine_data(self):
        """
        Function to extract meaningful results from compressed files
        """
        
        self.tot_pop = [[ ], [ ]]
        self.sim_time = [[ ], [ ]]
        self.ex = [[ ], [ ]]
        self.ey = [[ ], [ ]]
        self.ez = [[ ], [ ]]
        self.vesc = [[ ], [ ]]
        self.eKE = [[ ], [ ]]
        self.ePE = [[ ], [ ]]
        self.Nclose = [[ ], [ ]]
        self.ejecBool = [[ ], [ ]]

        ejec_data = natsort.natsorted(glob.glob('data/ejection_stats/*'))
        for file_ in range(len(ejec_data)):
            with open(ejec_data[file_], 'rb') as input_file:
                data_file = pkl.load(input_file)
                if data_file.iloc[0][0] == 'Hermite':
                    int_idx = 0
                else:
                    int_idx = 1

                self.tot_pop[int_idx].append(data_file.iloc[0][1])
                self.sim_time[int_idx].append(data_file.iloc[0][2])
                self.ex[int_idx].append(data_file.iloc[0][3])
                self.ey[int_idx].append(data_file.iloc[0][4])
                self.ez[int_idx].append(data_file.iloc[0][5])
                self.vesc[int_idx].append(data_file.iloc[0][6])
                self.eKE[int_idx].append(data_file.iloc[0][7])
                self.ePE[int_idx].append(data_file.iloc[0][8])
                self.Nclose[int_idx].append(data_file.iloc[0][9])
                self.ejecBool[int_idx].append(data_file.iloc[0][10])

        for int_ in range(2):
            self.ex[int_] = np.asarray(self.ex[int_], dtype = 'float')
            self.ey[int_] = np.asarray(self.ey[int_], dtype = 'float')
            self.ez[int_] = np.asarray(self.ez[int_], dtype = 'float')
            self.vesc[int_] = np.asarray(self.vesc[int_], dtype = 'float')
            self.eKE[int_] = np.asarray(self.eKE[int_], dtype = 'float')
            self.ePE[int_] = np.asarray(self.ePE[int_], dtype = 'float')

    def energy_plotters(self):
        """
        Function to plot KE vs. PE of the ejected particle
        """

        plot_ini = plotter_setup()
        MW_code = MWpotentialBovy2015()

        ejec_KEt = np.asarray(self.eKE)
        for int_ in range(1):
            ejec_KE = np.asarray(ejec_KEt[int_][np.isfinite(ejec_KEt[int_])])
            ejec_PE = np.asarray(self.ePE[int_][np.isfinite(ejec_KEt[int_])])
            ex = self.ex[int_][np.isfinite(self.ex[int_])]
            ey = self.ey[int_][np.isfinite(self.ey[int_])] 
            ez = self.ez[int_][np.isfinite(self.ez[int_])]

            ejected_KEb = np.asarray(ejec_KE[(ejec_KE < 1e46) & (abs(ejec_PE) < 1e46)])
            ejected_PEb = np.asarray(ejec_PE[(ejec_KE < 1e46) & (abs(ejec_PE) < 1e46)])
            final_posb = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(ex, ey, ez)])[(ejec_KE < 1e46) & (abs(ejec_PE) < 1e46)]
            final_posb[np.isnan(final_posb)] = 0
            ejected_KEb /= max(ejected_KEb)
            ejected_PEb /= max(ejected_KEb)

            PE_MW = ((1e3 | units.MSun) * MW_code.get_potential_at_point(0, 1 | units.pc, 1 | units.pc, 1 | units.pc)).value_in(units.J) / max(ejec_KE)

            ejected_KE = np.asarray([i/np.nanmax(ejec_KE) for i in ejec_KE])
            ejected_PE = np.asarray([i/np.nanmax(ejec_KE) for i in ejec_PE])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(ex, ey, ez)])
            final_pos[np.isnan(final_pos)] = 0

            linex = np.linspace(-0.01*np.min(ejected_KE), 10)
            liney = -linex

            plt.figure(figsize=(13, 4))
            gs = gridspec.GridSpec(1, 2)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0,1])
            for ax_ in [ax1, ax2]:
                ax_.plot(linex, liney, color = 'black', linestyle = '-.', zorder = 1)
                ax_.set_ylabel(r'$E_P/E_{{K, {\rm{max}}}}$')
                ax_.set_xlabel(r'$E_K/E_{{K, {\rm{max}}}}$')
                ax_.set_xlim(0, 1.1)
                ax_.set_ylim(-1.1, 0)
                plot_ini.tickers(ax_, 'plot')
                
            ax1.set_title(str(self.integrator[int_]))
            colour_axes = ax1.scatter(ejected_KE, ejected_PE, edgecolor = 'black', c=final_pos, s = 20, zorder = 2)
            cbar = plt.colorbar(colour_axes, ax=ax1, label = r'$r_{\rm{SMBH}}$ [pc]')

            bin2d_sim, xedge, yedge, img = ax2.hist2d(ejected_KE, ejected_PE, bins=(200, 200), range=([0,1.1],[1.1*min(ejected_PE),0]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [0, 1.1, 1.1*min(ejected_PE), 0]
            contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            plt.savefig('figures/ejection_stats/KEPE_diagram_'+str(self.integrator[int_])+'_.pdf', dpi=500, bbox_inches='tight')

    def vejec_plotter(self):
        """
        Function to plot the ejection velocity
        """

        plot_ini = plotter_setup()
        MW_code = MWpotentialBovy2015()
        vesc_MW = (np.sqrt(2)*MW_code.circular_velocity(0.1 | units.parsec) + np.sqrt(2*constants.G*(4e6 | units.MSun)/(0.1 | units.parsec))).value_in(units.kms)

        ymin = []
        ymax = []
        norm_min = np.log10(min(self.sim_time[0]))#, min(self.sim_time[1]))#, min(avg_surv[1])))
        norm_max = np.log10(max(self.sim_time[0]))#, max(self.sim_time[1]))
        normalise = plt.Normalize(norm_min, norm_max)

        fig = plt.figure(figsize=(15, 13))
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        
        ax = [ax1, ax2, ax3, ax4]
        for ax_ in [ax1]:#, ax3]:
            ax_.set_ylabel(r'$\langle v_{ej} \rangle$ [km/s]')
            ax_.axhline(vesc_MW, color = 'black', linestyle = ':')
            ax_.text(15, 660, r'$v_{\rm{esc, MW}}$')
            ax_.xaxis.labelpad = 30
            plot_ini.tickers(ax_, 'plot')
        for ax_ in [ax2]:#, ax4]:
            ax_.set_xlabel(r'$v_{ejec}$ [km s$^{-1}$]')
            ax_.set_ylabel(r'Frequency')
            plot_ini.tickers(ax_, 'plot')
            ax_.axvline(vesc_MW, linestyle = ':', color = 'black')
        ax_iter = -1

        ax_title = ['Ejection Velocity Histogram \n Hermite', '\n GRX']
        colours = ['red', 'blue']
        for int_ in range(1):
            ax_iter += 1
            tot_pop = np.asarray(self.tot_pop[int_])
            sim_time = np.asarray(self.sim_time[int_])
            vesc = np.asarray(self.vesc[int_])

            in_pop = np.unique(tot_pop[int_])
            avg_vesc = np.empty(len(in_pop))
            avg_surv = np.empty(len(in_pop))

            pops = [ ]
            samples = [ ]
            iter = -1
            for pop_ in in_pop:
                iter += 1
                pops.append(pop_)
                idx = np.where((tot_pop == pop_))[0]
                samples.append(len(vesc[idx]))
                avg_vesc[iter] = np.mean(vesc[idx])
                avg_surv[iter] = np.mean(sim_time[idx])

                ymin.append(min(avg_vesc))
                ymax.append(max(avg_vesc))

            colour_axes = ax[3*ax_iter].scatter(pops, avg_vesc, edgecolors='black', c = np.log10(avg_surv), norm = normalise, zorder = 3)
            ax[3*ax_iter].set_title(ax_title[int_])
            for j, xpos in enumerate(pops):
                ax[3*ax_iter].text(pops[j], -150, str(self.integrator[int_])+'\n'+str('{:.0f}'.format(samples[j])), fontsize = 'xx-small', ha = 'center' )
            plot_ini.tickers_pop(ax[ax_iter], pops)
            
            n, bins, patches = ax[3*ax_iter].hist(vesc, 20, histtype = 'step', color=colours[int_])
            n, bins, patches = ax[3*ax_iter].hist(vesc, 20, color=colours[int_], alpha = 0.4)

        """for ax_ in [ax1]:#, ax3]:
            ax_.set_ylim(min(ymin), max(ymax))
        for ax_ in [ax2]:#, ax4]:
            ax_.set_xlim(min(ymin), max(ymax))
            ax_.text(vesc_MW*(1+0.05), 0.9, r'$v_{esc, MW}$', rotation = 270, horizontalalignment = 'center')"""

        plt.savefig('figures/ejection_stats/vejection.pdf', dpi = 300, bbox_inches='tight')

class event_tracker(object):
    """
    Class to take stats of ejection vs. mergers
    """
    def __init__(self):

        plot_ini = plotter_setup()
        chaos_data = glob.glob('data/Hermite/no_addition/chaotic_simulation/*')
        chaos_data_GRX = glob.glob('data/GRX/no_addition/chaotic_simulation/*')
        chaos_data = [natsort.natsorted(chaos_data), natsort.natsorted(chaos_data_GRX)]

        init_pop = [[ ], [ ]]
        merger = [[ ], [ ]]
        in_pop = [[ ], [ ]]
        frac_merge = [[ ], [ ]]
        colours = ['red', 'blue']
        labels = ['Hermite', 'GRX']

        fig, ax = plt.subplots()
        ax.set_ylabel(r'$N_{\rm{merge}}/N_{\rm{sim}}$')
        ax.set_ylim(0,1.05)
        for int_ in range(2):
            for file_ in range(len(chaos_data[int_])):
                with open(chaos_data[int_][file_], 'rb') as input_file:
                    data = pkl.load(input_file)
                    init_pop[int_].append(len(data.iloc[0][8])+1)
                    merger[int_].append(data.iloc[0][10])

            in_pop[int_] = np.unique(init_pop[int_])
            iter = -1
            for pop_ in in_pop[int_]:
                iter +=1
                indices = np.where((init_pop[int_] == pop_))[0]
                temp_frac = [merger[int_][i] for i in indices]
                frac_merge[int_].append(np.mean(temp_frac))

            ax.scatter(in_pop[int_], frac_merge[int_], color = colours[int_], edgecolors = 'black', label = labels[int_])
        plot_ini.tickers_pop(ax, in_pop[0])
        plt.legend()
        plt.savefig('figures/ejection_stats/SMBH_merge_fraction.pdf', dpi=300, bbox_inches='tight')
