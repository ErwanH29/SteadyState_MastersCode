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
        self.chaos_ini_parti_data, self.chaos_fparti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, \
        self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_init_mass_data, \
        self.chaos_inj_mass_data, self.chaos_eje_mass_data = stats_chaos_extractor(dirC)  

        return self.chaos_ini_parti_data, self.chaos_fparti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, \
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
        pop, psamples = np.unique(filt_finparti, return_counts = True)

        return pop, psamples, filt_finparti, filt_stabtime

    def mean_plots(self, ax, ydata, pops):
        """
        Function to set up the various stability time plotters
        
        Inputs:
        ax:     Array containing axes in plots
        ydata:  The y-values for the plot
        """

        plot_ini = plotter_setup()
        mtick_formatter = [mtick.FormatStrFormatter('%0.1f'), mtick.FormatStrFormatter('%0.0f')]
        ylims = [[0, 1.2*(ydata)], [0, 105]]
    
        for i in range(len(ax)):
            plot_ini.tickers_pop(ax[i], pops)
            ax[i].yaxis.set_major_formatter(mtick_formatter[i])
            ax[i].set_ylim(ylims[i])

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

        chaos_ini_parti_data, chaos_fparti_data, chaos_number_mergers, chaos_cumulative_mm, chaos_simulated_end, \
        chaos_ejected_parti, chaos_stab_time_data, chaos_init_dist_data, chaos_init_mass_data, \
        chaos_inj_mass_data, chaos_eje_mass_data = self.chaos_extract(dirH)

        chaos_ini_parti_dataG, chaos_fparti_dataG, chaos_number_mergersG, chaos_cumulative_mmG, chaos_simulated_endG, \
        chaos_ejected_partiG, chaos_stab_time_dataG, chaos_init_dist_dataG, chaos_init_mass_dataG, \
        chaos_inj_mass_dataG, chaos_eje_mass_dataG = self.chaos_extract(dirG)

        in_dist = np.unique(chaos_init_dist_data)
        in_mass = np.unique(chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            init_dist_idx_chaos = np.where((np.asarray(chaos_init_dist_data) == dist_))
            init_dist_idx_chaosG = np.where((np.asarray(chaos_init_dist_dataG) == dist_))

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            tot_pop  = [ ]   
            iter = -1   
            coloursycler = cycle(colour_picker())
            for mass_ in in_mass:
                iter += 1
                pop, psamples, fparti, stab_time = self.index_extractor(chaos_init_mass_data, chaos_fparti_data, 
                                                                        chaos_stab_time_data, init_dist_idx_chaos, mass_)
                popG, psamplesG, fpartiG, stab_timeG = self.index_extractor(chaos_init_mass_dataG, chaos_fparti_dataG, 
                                                                            chaos_stab_time_dataG, init_dist_idx_chaosG, mass_)

                N_parti_avg = [ ]
                pop_id = np.argwhere(pop > 2)
                pop = pop[pop_id]
                psamples = psamples[pop_id]

                N_parti_avgG = [ ]
                pop_idG = np.argwhere(popG > 2)
                popG = popG[pop_idG]
                psamplesG = psamplesG[pop_idG]

                tot_pop.append(max(pop))

                full_simul = []
                for pop_, samp_ in zip(pop, psamples):
                    N_parti = np.argwhere(fparti == pop_)
                    N_parti_avg.append(np.mean(stab_time[N_parti]))
                    idx = np.where(stab_time[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_time[N_parti])
                    full_simul.append(ratio)
                N_parti_avg = np.asarray(N_parti_avg)
                print('For Hermite, # of full simulations per population: ', pop.flatten(), full_simul)

                full_simul = []
                for pop_, samp_ in zip(popG, psamplesG):
                    N_parti = np.argwhere(fpartiG == pop_)
                    N_parti_avgG.append(np.mean(stab_timeG[N_parti]))
                    idx = np.where(stab_timeG[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_timeG[N_parti])
                    full_simul.append(ratio)
                N_parti_avgG = np.asarray(N_parti_avgG)
                print('For GRX, # of full simulations per population: ', popG.flatten(), full_simul)
                
                ax.scatter(pop, N_parti_avg, color = 'red', edgecolor = 'black', zorder = 3,
                           label = r'Hermite')                
                ax.scatter(popG, N_parti_avgG, edgecolor = 'black', color = 'blue', 
                           zorder = 3, label = r'Hermite GRX')    
                ax.set_ylabel(r'$t_{surv}$ [Myr]')   

            for j, xpos in enumerate(pop):
                ax.text(pop[j][0], -0.1*max(N_parti_avg), '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(psamples[j][0])), fontsize = 'xx-small', ha = 'center' )
            for j, xpos in enumerate(popG):
                ax.text(popG[j][0], -0.11*max(N_parti_avg), 'GRX: '+str('{:.0f}'.format(psamplesG[j][0])), fontsize = 'xx-small', ha = 'center' )

            pop = np.array([float(i) for i in pop])
            N_parti_avg = np.array([ float(i) for i in N_parti_avg])   
             
            self.mean_plots([ax], max(N_parti_avg), pop)
            ax.set_xlim(5,105)

            p0 = (6,  0.50)
            params, cv = scipy.optimize.curve_fit(log_fit, pop, N_parti_avg, p0)
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
                fparti = chaos_fparti_data[init_dist_idx]
                stab_time = chaos_stab_time_data[init_dist_idx]

                tot_pop = [ ]
                N_parti_avg = [ ]
                std = [ ]
                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column

                fparti = fparti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop, psamples = np.unique(fparti, return_counts=True)
                pop_id = np.argwhere(pop > 2)
                pop = pop[pop_id]

                if len(pop) == 0:
                    pass
                else:
                    psamples = psamples[pop_id]
                    tot_pop.append(max(pop))

                    for pop_, samp_ in zip(pop, psamples):
                        N_parti = np.argwhere(fparti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))
                    y_max.append(max(np.add(N_parti_avg, std)))

                    ax.errorbar(pop, N_parti_avg, color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop, np.add(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop, np.subtract(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop, N_parti_avg, color = 'black')
                    
                    ax.text(85, 0.96*(max(np.add(N_parti_avg, std))), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} =$ '+str(mass_[0])+r' $M_\odot$')
                    ax.set_xlim(5, 105)
                    ax.set_ylim(0, 1.1*max(np.add(N_parti_avg, std)))
                    ax.set_ylabel(r'$t_{eject}$ [Myr]')
                    ax.set_title(r'Spread in Stability Time')
                    plot_ini.tickers(ax, 'plot')
                    plt.savefig('figures/steady_time/const_pop_chaotic_stab_time_equal_dist_'+str(dist_)+'_err_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')
