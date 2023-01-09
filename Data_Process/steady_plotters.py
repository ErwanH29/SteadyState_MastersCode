from amuse.lab import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import iqr
import scipy.optimize
from scipy.optimize import OptimizeWarning
from scipy.interpolate import make_interp_spline
import pickle as pkl

class stability_plotters(object):
    """
    Class to extract data and setup plots for simulations where IMBH particles are added over time
    """
        
    def extract(self, dirC):
        """
        Function to extract data from simulations who end in ejection
        """

        warnings.filterwarnings("ignore", category=OptimizeWarning) 

        self.fparti_data, self.stab_time_data, self.init_dist_data, self.imass_data = stats_chaos_extractor(dirC)  

        return self.fparti_data, self.stab_time_data, self.init_dist_data, self.imass_data

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
        pop, psamp = np.unique(filt_finparti, return_counts = True)

        return pop, psamp, filt_finparti, filt_stabtime
    
    def overall_steady_plotter(self):
        """
        Function to plot stability time for constant distances
        """
            
        def log_fit(xval, slope, beta, log_c, beta2, yint):
            return slope *(xval**0.5*np.log(log_c*xval))**beta
            
        plot_ini = plotter_setup()
        dirH = 'data/Hermite/no_addition/chaotic_simulation/*'
        dirG = 'data/GRX/no_addition/chaotic_simulation/*'

        fparti_data, stab_time_data, init_dist_data, imass_data = self.extract(dirH)
        fparti_dataG, stab_time_dataG, init_dist_dataG, imass_dataG = self.extract(dirG)

        idist_idx_chaos = np.where((np.asarray(init_dist_data) == init_dist_data))
        idist_idx_chaosG = np.where((np.asarray(init_dist_dataG) == init_dist_dataG))

        tot_pop  = [ ]   
        pop, psamp, fparti, stab_time = self.index_extractor(imass_data, fparti_data, 
                                                             stab_time_data, idist_idx_chaos, 
                                                             imass_data)
        popG, psampG, fpartiG, stab_timeG = self.index_extractor(imass_dataG, fparti_dataG, 
                                                                 stab_time_dataG, idist_idx_chaosG, 
                                                                 imass_dataG)
        
        psampG = psampG[popG > 5]
        popG = popG[popG > 5]

        N_parti_avg = [[ ], [ ]]
        N_parti_std = [[ ], [ ]]
        avg_deviate = [[ ], [ ]]
        full_simul = [[ ], [ ]]
        fparti = [fparti, fpartiG]
        stab_time = [stab_time, stab_timeG]
        pop_id = np.argwhere(pop > 2)
        pop_idG = np.argwhere(popG > 2)
        pop = [pop[pop_id], popG[pop_idG]]
        psamp = [psamp[pop_id], psampG[pop_idG]]

        tot_pop.append(max(max(pop[0]), max(pop[1])))

        colors = ['red', 'blue']
        integrator = ['Hermite', 'GRX']
        std_max = [[ ], [ ]]
        std_min = [[ ], [ ]]
        
        temp_data = []
        for int_ in range(2):
            for pop_, samp_ in zip(pop[int_][pop[int_] > 5], psamp[int_][pop[int_] > 5]):
                N_parti = np.argwhere(fparti[int_] == pop_)
                time_arr = stab_time[int_][N_parti]
                stability = np.mean(time_arr)

                if int_ == 0:
                    if pop_ == 60 or pop_ == 70:
                        temp_data.append(stab_time[int_][N_parti])

                N_parti_avg[int_].append(stability)
                N_parti_std[int_].append((np.std(time_arr)))
                q1, q3 = np.percentile(time_arr, [25, 75])
                std_max[int_].append(q3)
                std_min[int_].append(q1)

                idx = np.where(time_arr == 100)[0]
                ratio = len(idx)/len(time_arr)
                full_simul[int_].append(ratio)
                avg_deviate[int_].append(abs(np.mean(time_arr-np.std(time_arr))))

            N_parti_avg[int_] = np.asarray(N_parti_avg[int_])
            N_parti_std[int_] = np.asarray(N_parti_std[int_])
            std_max[int_] = np.asarray(std_max[int_])
            std_min[int_] = np.asarray(std_min[int_])
            avg_deviate[int_] = np.asarray(avg_deviate[int_])
            pop[int_] = np.array([float(i) for i in pop[int_]])
            
        fig, ax = plt.subplots()
        hist_tails = np.concatenate((temp_data[0], temp_data[1]))
        hist_temp = [N_parti_avg[0][5], N_parti_avg[0][6]]
        ax.set_title(r'$\langle t_{\rm{surv}} \rangle$ for $N = 60, N = 70$')
        ax.axvline((4.87), color = 'black', linestyle = ':')
        ax.text((4.95), 17, r'$\langle t_{\rm{surv}}\rangle$', rotation = 270)
        ax.set_xlabel(r'Time [Myr]')
        ax.set_ylabel(r'Counts')
        n1, bins, patches = ax.hist((hist_tails), 30, histtype='step', color = 'black')
        n1, bins, patches = ax.hist((hist_tails), 30, alpha = 0.3, color = 'black')
        plot_ini.tickers(ax, 'plot')
        plt.savefig('figures/steady_time/stab_time_hist_6070.pdf', dpi = 300, bbox_inches='tight')
        plt.clf()

        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel(r'$\log_{10} t_{\rm{surv}}$ [Myr]') 
        plot_ini.tickers_pop(ax1, pop[0], 'Hermite')
        ax1.set_xlim(5,105)
        ax1.set_ylim(-1.5, 2.3)
        ax1.axhline(np.log10(259))
        ax1.axhline(np.log10(882))

        for int_ in range(2):
            N_parti_avg[int_] = np.array([float(i) for i in N_parti_avg[int_]])
            for j, xpos in enumerate(pop[int_][pop[int_] % 10 == 0]):
                if j == 0:
                    ax1.scatter(pop[int_][pop[int_] % 10 == 0], np.log10(N_parti_avg[int_][pop[int_] % 10 == 0]), color = colors[int_], edgecolor = 'black', zorder = 2, label = integrator[int_])
                else:
                    ax1.scatter(pop[int_][pop[int_] % 10 == 0], np.log10(N_parti_avg[int_][pop[int_] % 10 == 0]), color = colors[int_], edgecolor = 'black', zorder = 2)
            ax1.scatter(pop[int_][pop[int_] % 10 == 0], np.log10(std_min[int_][pop[int_] % 10 == 0]), color = colors[int_], marker = '_')
            ax1.scatter(pop[int_][pop[int_] % 10 == 0], np.log10(std_max[int_][pop[int_] % 10 == 0]), color = colors[int_], marker = '_')
            ax1.plot([pop[int_][pop[int_] % 10 == 0], pop[int_][pop[int_] % 10 == 0]], [np.log10(std_min[int_][pop[int_] % 10 == 0]), np.log10(std_max[int_][pop[int_] % 10 == 0])], color = colors[int_], zorder = 1)

        p0 = (100, -5, 20, 0.5, 60)
        params, cv = scipy.optimize.curve_fit(log_fit, pop[1][pop[1] % 10 == 0], (N_parti_avg[1][pop[1] % 10 == 0]), p0, maxfev = 10000, method = 'trf')
        slope, beta, log_c, beta2, y = params
        
        xtemp = np.linspace(10, 100, 1000)
        curve = [(log_fit(i, slope, beta, log_c, beta2, y)) for i in xtemp]

        ax1.plot(xtemp, np.log10(curve), zorder = 1, color = 'black', ls = '-.')
        ax1.legend()
        plt.savefig('figures/steady_time/stab_time_mean.pdf', dpi = 300, bbox_inches='tight')
        plt.clf()

        print('Factor:       ', slope)
        print('log Factor:   ', log_c)
        print('Power Factor: ', beta)
      
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel(r'$\log_{10} t_{\rm{surv}}$ [Myr]') 
        ax1.set_xlim(8,55)
        ax1.set_ylim(-1.3, 2.3)
        xints = [i for i in range(8, 1+int(max(pop[1]))) if i % 5 == 0]
        ax1.set_xlabel(r'IMBH Population [$N$]')
        ax1.yaxis.set_ticks_position('both')
        ax1.xaxis.set_ticks_position('both')
        ax1.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax1.tick_params(axis="y", which = 'both', direction="in")
        ax1.tick_params(axis="x", which = 'both', direction="in")     
        ax1.set_xticks(xints)
        
        N_parti_avg[int_] = np.array([float(i) for i in N_parti_avg[int_]])
        xtemp = np.linspace(10, 35, 1000)
        curve = [(log_fit(i, slope, beta, log_c, beta2, y)) for i in xtemp]
        for j, xpos in enumerate(pop[1]):
            if pop[1][j] > 5 and pop[1][j] <= 40:
                if j == 0:
                    ax1.scatter(pop[1][pop[1] < 40], np.log10(N_parti_avg[1][pop[1] < 40]), color = colors[1], edgecolor = 'black', zorder = 2, label = integrator[1])
                else:
                    ax1.scatter(pop[1][pop[1] < 40], np.log10(N_parti_avg[1][pop[1] < 40]), color = colors[1], edgecolor = 'black', zorder = 2)
        ax1.scatter(pop[1][pop[1] < 40], np.log10(std_min[1][pop[1] < 40]), color = colors[1], marker = '_')
        ax1.scatter(pop[1][pop[1] < 40], np.log10(std_max[1][pop[1] < 40]), color = colors[1], marker = '_')
        ax1.plot([pop[1][pop[1] < 40], pop[1][pop[1] < 40]], [np.log10(std_min[1][pop[1] < 40]), np.log10(std_max[1][pop[1] < 40])], color = colors[1], zorder = 1)
        print(pop, psamp)

        ax1.plot(xtemp, np.log10(curve), zorder = 1, color = 'black', ls = '-.')
        plot_ini.tickers_pop(ax1, pop[0], 'GRX')
        plt.savefig('figures/steady_time/stab_time_mean_GRX.pdf', dpi = 300, bbox_inches='tight')

        fig = plt.figure(figsize=(15, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax = [ax1, ax2]
        for int_ in range(2):
            x_arr = np.linspace(10, max(pop[int_]), 100)
            smooth_curve = make_interp_spline(pop[int_], avg_deviate[int_])
            ax[int_].set_ylabel(r'$\langle (t_{\rm{surv}} - \sigma_{\rm{surv}}) \rangle$')
            ax[int_].set_title(integrator[int_])
            ax[int_].plot(x_arr, smooth_curve(x_arr), color = colors[int_], zorder = 1)
            ax[int_].scatter(pop[int_][pop[int_] > 5], avg_deviate[int_][pop[int_] > 5], color = colors[int_], edgecolors='black', zorder = 2)   
            ax[int_].set_ylim(0,1.05*max(avg_deviate[int_][pop[int_] > 5]))
        plot_ini.tickers_pop(ax1, pop[0], 'Hermite')
        plot_ini.tickers_pop(ax2, pop[1], 'GRX')
        ax1.set_xlim(6,105)
        ax2.set_xlim(6,55)
        plt.savefig('figures/steady_time/stab_time_residuals.pdf', dpi = 300, bbox_inches='tight')
        plt.clf()

        with open('figures/steady_time/Sim_summary.txt', 'w') as file:
            for int_ in range(2):
                file.write('\n\nFor'+str(integrator[int_])+', # of full simulations per population:  '+str(pop[int_].flatten()))
                file.write('\n                                              '+str(full_simul[int_]))
                file.write('\nThe slope of the curve goes as:               '+str(slope))
                file.write('\nThe power-law of the lnN goes as:             '+str(beta))
                file.write('\nThe logarithmic factor goes as:               '+str(log_c))
                file.write('\nThe final raw data:                           '+str(pop[int_].flatten()))
                file.write('\nSimulated time [Myr]                          '+str(N_parti_avg[int_].flatten()))

