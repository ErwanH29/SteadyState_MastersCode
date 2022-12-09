from amuse.lab import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import iqr
import scipy.optimize
from scipy.optimize import OptimizeWarning
from scipy.interpolate import make_interp_spline

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

        def f(x, mu, sigma) :
            return 1/(np.sqrt(2*np.pi)*sigma*x)*np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))

        def log_slope(slope):
            rclust = 0.3 | units.pc
            Mclust = 10**4 | units.MSun
            slope = slope * (1 | units.Myr)
            trlx_coeff = (rclust/(constants.G*Mclust)**1/3)

            return np.log(abs(slope.value_in(units.s)))/np.log(trlx_coeff.value_in(units.s**2/units.m**2))
            
        def log_fit(xval, slope, beta, log_c, y_int):
            return slope * (1/np.log(log_c*xval))**beta# + y_int
            
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
        for int_ in range(2):
            for pop_, samp_ in zip(pop[int_], psamp[int_]):
                N_parti = np.argwhere(fparti[int_] == pop_)
                time_arr = stab_time[int_][N_parti]
                stability = np.mean(time_arr)

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
            pop[int_] = np.array([float(i) for i in pop[int_]])
        
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel(r'$\log_{10} t_{\rm{surv}}$ [Myr]') 
        plot_ini.tickers_pop(ax1, pop[0])
        ax1.set_xlim(5,105)
        ax1.xaxis.labelpad = 30
        ax1.set_ylim(-2.5, 2.3)

        for int_ in range(2):
            N_parti_avg[int_] = np.array([float(i) for i in N_parti_avg[int_]])
            for j, xpos in enumerate(pop[int_]):
                ax1.text(pop[int_][j], -2.75*(1+0.04*int_), str(integrator[int_])+': '+str('{:.0f}'.format(psamp[int_][j][0])), fontsize = 'xx-small', ha = 'center' )
                if j == 0:
                    ax1.scatter(pop[int_], np.log10(N_parti_avg[int_]), color = colors[int_], edgecolor = 'black', zorder = 2, label = integrator[int_])
                else:
                    ax1.scatter(pop[int_], np.log10(N_parti_avg[int_]), color = colors[int_], edgecolor = 'black', zorder = 2)
            ax1.scatter(pop[int_], np.log10(std_min[int_]), color = colors[int_], marker = '_')
            ax1.scatter(pop[int_], np.log10(std_max[int_]), color = colors[int_], marker = '_')
            ax1.plot([pop[int_], pop[int_]], [np.log10(std_min[int_]), np.log10(std_max[int_])], color = colors[int_], zorder = 1)

        p0 = (100, -5, 20, 40)
        params, cv = scipy.optimize.curve_fit(log_fit, pop[1][:-1], (N_parti_avg[1][:-1]), p0, maxfev = 10000, method = 'trf')
        slope, beta, log_c, yint = params
        
        slope_str = str('{:.2f}'.format(slope))
        logc_str = str('{:.2f}'.format(log_c))
        beta_str = str('{:.2f}'.format(beta))
        xtemp = np.linspace(10, 100, 1000)
        curve = [(log_fit(i, slope, beta, log_c, yint)) for i in xtemp]

        ax1.plot(xtemp, np.log10(curve), zorder = 1, color = 'black', ls = '-.')
        ax1.text(72, 1.5, r'$t_{{\rm surv}} \approx{{{}}}(\frac{{N}}{{\ln({{{}N}})}}$'.format(slope_str[:3], logc_str)+r'$)^{{{}}}$'.format(beta_str)+' Myr')
        ax1.legend()
        plt.savefig('figures/steady_time/stab_time_mean.pdf', dpi = 300, bbox_inches='tight')

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
            ax[int_].scatter(pop[int_], avg_deviate[int_], color = colors[int_], edgecolors='black', zorder = 2)   
            plot_ini.tickers_pop(ax[int_], pop[0])
        plt.savefig('figures/steady_time/stab_time_residuals.pdf', dpi = 300, bbox_inches='tight')
        plt.clf()

        with open('figures/steady_time/Sim_summary.txt', 'w') as file:
            for int_ in range(2):
                file.write('\n\nFor'+str(integrator[int_])+', # of full simulations per population:  '+str(pop[int_].flatten()))
                file.write('\n                                                   '+str(full_simul[int_]))
                file.write('\nThe slope of the full curve goes as:               '+str(log_slope(slope)))
                #file.write('\nwith errors:                                       '+str(errs_upper)+' '+str(errs_lower))
                file.write('\nThe logarithmic prefactor goes as:                 '+str(log_c))
                file.write('\nThe final raw data:                                '+str(pop[int_].flatten()))
                file.write('\nSimulated time [Myr]                               '+str(N_parti_avg[int_].flatten()))
                file.write('\nStandard dev. [Myr]:                               '+str(N_parti_std[int_].flatten()))



print('...steady_plotter...')
cst = stability_plotters()
cst.overall_steady_plotter()