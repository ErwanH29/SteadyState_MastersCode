from amuse.lab import *
from file_logistics import *
from spatial_plotters import *
from itertools import cycle
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy.optimize

class stability_plotters(object):
    """
    Class to extract data and setup plots for simulations where IMBH particles are added over time
    """
        
    def chaos_extract(self, dirC):
        """
        Function to extract data from simulations who end in ejection
        """
        self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, \
        self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_init_mass_data, \
        self.chaos_inj_mass_data, self.chaos_eje_mass_data = stats_chaos_extractor(dirC)  

        return self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, \
               self.chaos_simulated_end, self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, \
               self.chaos_init_mass_data, self.chaos_inj_mass_data, self.chaos_eje_mass_data

    def index_extractor(self, init_mass, final_parti_data, stab_time_data, indices, vals):
        """
        Function to extract indices of depending on the wanted plot (function of mass/distance)
        
        Inputs:
        init_mass:         The set of initial masses
        final_parti_data:  The data of the final configuration
        stab_time_data:    The stability time
        indices:           The indices corresponding to the fashion you wish to separate the plots
        vals:              The variable which will compose different coloured lines in the plot
        """
        
        mass_array = init_mass[indices]
        final_part = final_parti_data[indices]
        stab_time = stab_time_data[indices]

        filtered_m = np.where((mass_array == vals).all(1))[0]
        filt_finparti = final_part[filtered_m]
        filt_stabtime = stab_time[filtered_m]
        pop_size, pop_samples = np.unique(filt_finparti, return_counts = True)

        return pop_size, pop_samples, filt_finparti, filt_stabtime

    def mean_plots(self, axis, ydata, pops):
        """
        Function to set up the various stability time plotters
        
        Inputs:
        ax:     Array containing axes in plots
        ydata:  The y-values for the plot
        """

        plot_ini = plotter_setup()
        mtick_formatter = [mtick.FormatStrFormatter('%0.1f'), mtick.FormatStrFormatter('%0.0f')]
        ylims = [[0, 1.2*(ydata)], [0, 105]]
    
        for i in range(len(axis)):
            plot_ini.tickers_pop(axis[i], pops)
            axis[i].yaxis.set_major_formatter(mtick_formatter[i])
            axis[i].set_ylim(ylims[i])

    def stable_extract(self, dirS):
        """
        Function to extract data  from simulations who end with stable state
        """
        self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
        self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data = stats_stable_extractor(dirS)

        return self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
               self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data

    def massdep_plotter(self):
        """
        Function to plot stability time for constant distances
        """

        def log_fit(xval, slope, yint):
            return (slope) /( (xval)*np.log(xval)) + yint

        plot_ini = plotter_setup()

        dirH = 'data/Hermite/no_addition/chaotic_simulation/*'
        dirG = 'data/GRX/no_addition/chaotic_simulation/*'

        chaos_ini_parti_data, chaos_fin_parti_data, chaos_number_mergers, chaos_cumulative_mm, chaos_simulated_end, \
        chaos_ejected_parti, chaos_stab_time_data, chaos_init_dist_data, chaos_init_mass_data, \
        chaos_inj_mass_data, chaos_eje_mass_data = self.chaos_extract(dirH)

        chaos_ini_parti_data_GRX, chaos_fin_parti_data_GRX, chaos_number_mergers_GRX, chaos_cumulative_mm_GRX, chaos_simulated_end_GRX, \
        chaos_ejected_parti_GRX, chaos_stab_time_data_GRX, chaos_init_dist_data_GRX, chaos_init_mass_data_GRX, \
        chaos_inj_mass_data_GRX, chaos_eje_mass_data_GRX = self.chaos_extract(dirG)

        in_dist = np.unique(chaos_init_dist_data)
        in_mass = np.unique(chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            init_dist_idx_chaos = np.where((np.asarray(chaos_init_dist_data) == dist_))
            init_dist_idx_chaos_GRX = np.where((np.asarray(chaos_init_dist_data_GRX) == dist_))

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            tot_pop  = [ ]   
            iter = -1   
            coloursycler = cycle(colour_picker())
            for mass_ in in_mass:
                iter += 1
                pop_size, pop_samples, fin_parti, stab_time = self.index_extractor(chaos_init_mass_data, chaos_fin_parti_data, 
                                                                                   chaos_stab_time_data, init_dist_idx_chaos, mass_)
                pop_size_GRX, pop_samples_GRX, fin_parti_GRX, stab_time_GRX = self.index_extractor(chaos_init_mass_data_GRX, chaos_fin_parti_data_GRX, 
                                                                                                   chaos_stab_time_data_GRX, init_dist_idx_chaos_GRX, mass_)

                N_parti_avg = [ ]
                pop_id = np.argwhere(pop_size > 2)
                pop_size = pop_size[pop_id]
                pop_samples = pop_samples[pop_id]

                N_parti_avg_GRX = [ ]
                pop_id_GRX = np.argwhere(pop_size_GRX > 2)
                pop_size_GRX = pop_size_GRX[pop_id_GRX]
                pop_samples_GRX = pop_samples_GRX[pop_id_GRX]

                tot_pop.append(max(pop_size))

                full_simul = []
                for pop_, samp_ in zip(pop_size, pop_samples):
                    N_parti = np.argwhere(fin_parti == pop_)
                    N_parti_avg.append(np.mean(stab_time[N_parti]))
                    idx = np.where(stab_time[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_time[N_parti])
                    full_simul.append(ratio)
                N_parti_avg = np.asarray(N_parti_avg)
                print('For Hermite, # of full simulations per population: ', pop_size.flatten(), full_simul)

                full_simul = []
                for pop_, samp_ in zip(pop_size_GRX, pop_samples_GRX):
                    N_parti = np.argwhere(fin_parti_GRX == pop_)
                    N_parti_avg_GRX.append(np.mean(stab_time_GRX[N_parti]))
                    idx = np.where(stab_time_GRX[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_time_GRX[N_parti])
                    full_simul.append(ratio)
                N_parti_avg_GRX = np.asarray(N_parti_avg_GRX)
                print('For GRX, # of full simulations per population: ', pop_size_GRX.flatten(), full_simul)
                
                ax.scatter(pop_size, N_parti_avg, color = 'red', edgecolor = 'black', zorder = 3,
                           label = r'Hermite')                
                ax.scatter(pop_size_GRX, N_parti_avg_GRX, edgecolor = 'black', color = 'blue', 
                           zorder = 3, label = r'Hermite GRX')    
                ax.set_ylabel(r'$t_{surv}$ [Myr]')   

            for j, xpos in enumerate(pop_size):
                ax.text(pop_size[j][0], -0.1*max(N_parti_avg), '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
            for j, xpos in enumerate(pop_size_GRX):
                ax.text(pop_size_GRX[j][0], -0.1*max(N_parti_avg_GRX), '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(pop_samples_GRX[j][0])), fontsize = 'xx-small', ha = 'center' )

            pop_size = np.array([float(i) for i in pop_size])
            N_parti_avg = np.array([ float(i) for i in N_parti_avg])   
             
            self.mean_plots([ax], max(N_parti_avg), pop_size)
            ax.set_xlim(5,105)

            p0 = (6,  0.50)
            params, cv = scipy.optimize.curve_fit(log_fit, pop_size, N_parti_avg, p0)
            slope, intercept = params
            red_slope = str('{:.2f}'.format(slope))
            xtemp = np.linspace(8, 105, 1000)
            ytemp = [log_fit(i, slope, intercept) for i in xtemp]
            
            ax.plot(xtemp, ytemp, zorder = 1, color = 'black', ls = '-.')
            ax.text(8, 4, r'$t_{{surv}} \approx \frac{{{}}}{{N\lnN}}$'.format(red_slope)+ ' Myr')
            ax.legend()
            ax.xaxis.labelpad = 30
            plt.savefig('figures/steady_time/const_population_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for dist_ in in_dist:
            for mass_ in in_mass: 
                y_max = []

                plt.clf()
                fig, ax = plt.subplots(figsize=(8,6))

                init_dist_idx = np.where((chaos_init_dist_data == dist_))
                init_mass = chaos_init_mass_data[init_dist_idx]
                fin_parti = chaos_fin_parti_data[init_dist_idx]
                stab_time = chaos_stab_time_data[init_dist_idx]

                tot_pop = [ ]
                N_parti_avg = [ ]
                std = [ ]
                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column

                fin_parti = fin_parti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)
                pop_id = np.argwhere(pop_size > 2)
                pop_size = pop_size[pop_id]

                if len(pop_size) == 0:
                    pass
                else:
                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))
                    y_max.append(max(np.add(N_parti_avg, std)))

                    ax.errorbar(pop_size, N_parti_avg, color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, N_parti_avg, color = 'black')
                    
                    ax.text(85, 0.96*(max(np.add(N_parti_avg, std))), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} =$ '+str(mass_[0])+r' $M_\odot$')
                    ax.set_xlim(5, 105)
                    ax.set_ylim(0, 1.1*max(np.add(N_parti_avg, std)))
                    ax.set_ylabel(r'$t_{eject}$ [Myr]')
                    ax.set_title(r'Spread in Stability Time')
                    plot_ini.tickers(ax, 'plot')
                    plt.savefig('figures/steady_time/const_pop_chaotic_stab_time_equal_dist_'+str(dist_)+'_err_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')
