from amuse.lab import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.optimize import OptimizeWarning

class stability_plotters(object):
    """
    Class to extract data and setup plots for simulations where IMBH particles are added over time
    """
        
    def extract(self, dirC):
        """
        Function to extract data from simulations who end in ejection
        """

        warnings.filterwarnings("ignore", category=OptimizeWarning) 

        self.iparti_data, self.fparti_data, self.no_mergers, self.cumuMM, self.sim_end, \
        self.ejec_parti, self.stab_time_data, self.init_dist_data, self.imass_data, \
        self.inj_mass_data, self.eje_mass_data = stats_chaos_extractor(dirC)  

        return self.iparti_data, self.fparti_data, self.no_mergers, self.cumuMM, \
               self.sim_end, self.ejec_parti, self.stab_time_data, self.init_dist_data, \
               self.imass_data, self.inj_mass_data, self.eje_mass_data

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

    def stable_extract(self, dirS):
        """
        Function to extract data  from simulations who end with stable state
        """
        self.iparti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
        self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data = stats_stable_extractor(dirS)

        return self.iparti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
               self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data
    
    def overall_steady_plotter(self):
        """
        Function to plot stability time for constant distances
        """

        def log_slope(slope):
            rclust = 0.3 | units.pc
            Mclust = 10**4 | units.MSun
            slope = slope * (1 | units.Myr)
            trlx_coeff = (rclust/(constants.G*Mclust)**1/3)

            return np.log(abs(slope.value_in(units.s)))/np.log(trlx_coeff.value_in(units.s**2/units.m**2))
            
        def log_fit(xval, slope, beta, log_c):
            return slope * ((xval)/np.log(log_c*xval))**beta
            
        plot_ini = plotter_setup()

        dirH = 'data/Hermite/no_addition/chaotic_simulation/*'
        dirG = 'data/GRX/no_addition/chaotic_simulation/*'

        iparti_data, fparti_data, no_mergers, cumuMM, sim_end, \
        ejec_parti, stab_time_data, init_dist_data, imass_data, \
        inj_mass_data, eje_mass_data = self.extract(dirH)

        iparti_dataG, fparti_dataG, no_mergersG, cumuMMG, sim_endG, \
        ejec_partiG, stab_time_dataG, init_dist_dataG, imass_dataG, \
        inj_mass_dataG, eje_mass_dataG = self.extract(dirG)

        idist_idx_chaos = np.where((np.asarray(init_dist_data) == init_dist_data))
        idist_idx_chaosG = np.where((np.asarray(init_dist_dataG) == init_dist_dataG))

        tot_pop  = [ ]   
        pop, psamp, fparti, stab_time = self.index_extractor(imass_data, fparti_data, 
                                                             stab_time_data, idist_idx_chaos, 
                                                             imass_data)

        popG, psampG, fpartiG, stab_timeG = self.index_extractor(imass_dataG, fparti_dataG, 
                                                                 stab_time_dataG, idist_idx_chaosG, 
                                                                 imass_dataG)

        N_parti_avg = [ ]
        N_parti_std = [ ]
        pop_id = np.argwhere(pop > 2)
        pop = pop[pop_id]
        psamp = psamp[pop_id]

        N_parti_avgG = [ ]
        N_parti_stdG = [ ]
        pop_idG = np.argwhere(popG > 2)
        popG = popG[pop_idG]
        psampG = psampG[pop_idG]

        tot_pop.append(max(pop))

        full_simul_H = []
        avg_deviate_H = []
        for pop_, samp_ in zip(pop, psamp):
            N_parti = np.argwhere(fparti == pop_)
            N_parti_avg.append(np.mean(stab_time[N_parti]))
            N_parti_std.append((np.std(stab_time[N_parti])))
            idx = np.where(stab_time[N_parti] == 100)[0]
            ratio = len(idx)/len(stab_time[N_parti])
            full_simul_H.append(ratio)
            avg_deviate_H.append(abs(np.mean(stab_time[N_parti]-np.std(stab_time[N_parti]))))
        N_parti_avg = np.asarray(N_parti_avg)
        N_parti_std = np.asarray(N_parti_std)

        full_simul_G = []
        avg_deviate_G = []
        for pop_, samp_ in zip(popG, psampG):
            N_parti = np.argwhere(fpartiG == pop_)
            N_parti_avgG.append(np.mean(stab_timeG[N_parti]))
            N_parti_stdG.append(np.std(stab_timeG[N_parti]))
            idx = np.where(stab_timeG[N_parti] == 100)[0]
            ratio = len(idx)/len(stab_timeG[N_parti])
            full_simul_G.append(ratio)
            avg_deviate_G.append(abs(np.mean(stab_timeG[N_parti]-np.std(stab_timeG[N_parti]))))
        N_parti_avgG = np.asarray(N_parti_avgG)
        N_parti_stdG = np.asarray(N_parti_stdG)

        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax1.set_ylabel(r'$\log_{10} t_{\rm{surv}}$ [Myr]') 
        plot_ini.tickers_pop(ax1, pop)
        ax1.set_xlim(5,105)
        ax1.xaxis.labelpad = 30
        ax1.set_ylim(-2.5, 2.3)

        for j, xpos in enumerate(pop):
            ax1.text(pop[j][0], -2.75, 'Hermite: '+str('{:.0f}'.format(psamp[j][0])), fontsize = 'xx-small', ha = 'center' )
        for j, xpos in enumerate(popG):
            ax1.text(popG[j][0], -2.9, 'GRX: '+str('{:.0f}'.format(psampG[j][0])), fontsize = 'xx-small', ha = 'center' )

        pop = np.array([float(i) for i in pop])
        N_parti_avg = np.array([float(i) for i in N_parti_avg])

        p0 = (100, 5, 10)
        params, cv = scipy.optimize.curve_fit(log_fit, pop, (N_parti_avg), p0, maxfev = 10000, method = 'trf')
        slope, beta, log_c = params
        
        slope_str = str('{:.0f}'.format(slope))
        logc_str = str('{:.2f}'.format(log_c))
        beta_str = str('{:.2f}'.format(beta))
        xtemp = np.linspace(10, 100, 1000)
        curve = [(log_fit(i, slope, beta, log_c)) for i in xtemp]

        iter = 0
        for pop_ in pop:
            y_val = np.log10(N_parti_avg[iter]+N_parti_std[iter])
            y_val = [-y_val, y_val]
            if iter == 0:
                ax1.scatter(pop_, np.log10(N_parti_avg[iter]), color = 'red', edgecolor = 'black', zorder = 2, label = 'Hermite')
            else:
                ax1.scatter(pop_, np.log10(N_parti_avg[iter]), color = 'red', edgecolor = 'black', zorder = 2)
            ax1.scatter(pop_, y_val[-1], color = 'red', marker = '_')
            ax1.scatter(pop_, y_val[0], color = 'red', marker = '_')
            ax1.plot([pop_, pop_], y_val, color = 'red', zorder = 1)
            iter += 1
        
        iter = 0
        for pop_ in popG:
            y_val = np.log10(N_parti_avgG[iter]+N_parti_stdG[iter])
            y_val = [-y_val, y_val]
            if iter == 0:
                ax1.scatter(pop_, np.log10(N_parti_avgG[iter]), color = 'blue', edgecolor = 'black', zorder = 2, label = 'GRX')
            else:
                ax1.scatter(pop_, np.log10(N_parti_avgG[iter]), color = 'blue', edgecolor = 'black', zorder = 2)
            ax1.scatter(pop_, y_val[-1], color = 'blue', marker = '_')
            ax1.scatter(pop_, y_val[0], color = 'blue', marker = '_')
            ax1.plot([pop_, pop_], y_val, color = 'blue', zorder = 1)
            iter += 1

        ax1.plot(xtemp, np.log10(curve), zorder = 1, color = 'black', ls = '-.')
        ax1.text(72, 1.5, r'$t_{{\rm surv}} \approx{{{}}}(\frac{{N}}{{\ln({{{}N}})}}$'.format(slope_str[:3], logc_str)+r'$)^{{{}}}$'.format(beta_str)+' Gyr')
        ax1.legend()
        plt.savefig('figures/steady_time/stab_time_mean.pdf', dpi = 300, bbox_inches='tight')

        slope_avg = np.log10(N_parti_avg[-1]/N_parti_avg[0]) / 90
        slope_min = np.log10((N_parti_avg[-1]+N_parti_std[-1])/(N_parti_avg[0]-N_parti_std[0])) / 90
        slope_max = np.log10(abs((N_parti_avg[-1]-N_parti_std[-1])/(N_parti_avg[0]+N_parti_std[0]))) / 90
        ratio_min = slope_min / slope_avg
        ratio_max = slope_max / slope_avg

        best_fitH = np.poly1d(np.polyfit(pop, avg_deviate_H, 5))
        x_arr = np.linspace(10, 100, 100)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        ax.plot(x_arr, best_fitH(x_arr), color = 'red', zorder = 1)
        ax.scatter(pop, avg_deviate_H, color = 'red', edgecolors='black', label = 'Hermite', zorder = 2)
        ax.scatter(popG, avg_deviate_G, color = 'blue', edgecolors='black', label = 'GRX', zorder = 2)
        ax.legend()
        ax.set_ylabel(r'$\langle (t_{\rm{surv}} - \sigma_{\rm{surv}}) \rangle$')   
        plot_ini.tickers_pop(ax, pop)
        plt.savefig('figures/steady_time/stab_time_residuals.pdf', dpi = 300, bbox_inches='tight')
        plt.clf()

        errs_upper = log_slope(ratio_max * slope) - log_slope(slope)
        errs_lower = log_slope(ratio_min * slope) - log_slope(slope)
        #erra_upper = log_cmax - log_c
        #erra_lower = log_cmin - log_c
        with open('figures/steady_time/Sim_summary.txt', 'w') as file:
            file.write('For Hermite, # of full simulations per population: '+str(pop.flatten()))
            file.write('\n                                                   '+str(full_simul_H))
            file.write('\nThe slope of the full curve goes as:               '+str(log_slope(slope)))
            file.write('\nwith errors:                                       '+str(errs_upper)+' '+str(errs_lower))
            file.write('\n\nThe logarithmic prefactor goes as:                 '+str(log_c))
            #file.write('\nwith errors:                                       '+str(erra_upper)+' '+str(erra_lower))
            file.write('\n\n\nFor GRX, # of full simulations per population: '+str(pop.flatten()))
            file.write('\n                                               '+str(full_simul_G))
            #file.write('\nand the slope of the curve goes as:                '+str(np.log(slopeG)*np.exp(1)/10))
