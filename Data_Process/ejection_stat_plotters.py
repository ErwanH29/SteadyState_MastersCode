import matplotlib.pyplot as plt
from file_logistics import *
from spatial_plotters import *
from numpy.polynomial.polynomial import polyfit
import matplotlib.gridspec as gridspec

class data_setup(object):
    def __init__(self):
        """
        Extracts the required data
        """

        self.IMBH_tracker = bulk_stat_extractor('data/Hermite/particle_trajectory/*', 'Y')
        self.ejec_data = bulk_stat_extractor('data/Hermite/no_addition/chaotic_simulation/*', 'N')
        
        self.ex = np.empty((len(self.ejec_data)))
        self.ey = np.empty((len(self.ejec_data)))
        self.ez = np.empty((len(self.ejec_data)))
        self.vesc = np.empty((len(self.ejec_data)))

        self.ejec_KE = np.empty((len(self.ejec_data)))
        self.ejec_PE = np.empty((len(self.ejec_data)))
        self.Nclose = np.empty((len(self.ejec_data)))

        self.tot_mass = np.empty((len(self.ejec_data)))
        self.tot_pop = np.empty((len(self.ejec_data)))
        self.surv_time = np.empty((len(self.ejec_data)))
        self.ejected = np.empty(len(self.ejec_data))

        self.IMBH_tracker_GRX = bulk_stat_extractor('data/GRX/particle_trajectory/*', 'Y')
        self.ejec_data_GRX = bulk_stat_extractor('data/GRX/no_addition/chaotic_simulation/*', 'N')

        self.ex_GRX = np.empty((len(self.ejec_data_GRX)))
        self.ey_GRX = np.empty((len(self.ejec_data_GRX)))
        self.ez_GRX = np.empty((len(self.ejec_data_GRX)))
        self.vesc_GRX = np.empty((len(self.ejec_data_GRX)))

        self.ejec_KE_GRX = np.empty((len(self.ejec_data_GRX)))
        self.ejec_PE_GRX = np.empty((len(self.ejec_data_GRX)))
        self.Nclose_GRX = np.empty((len(self.ejec_data_GRX)))

        self.tot_mass_GRX = np.empty((len(self.ejec_data_GRX)))
        self.tot_pop_GRX = np.empty((len(self.ejec_data_GRX)))
        self.surv_time_GRX = np.empty((len(self.ejec_data_GRX)))
        self.ejected_GRX = np.empty(len(self.ejec_data_GRX))

class KE_PE_plotters(object):
    """
    Class to plot KE vs. PE
    """

    def energy_data_extract(self, filter):
        """
        Function to extract particle energy data
        
        Inputs:
        filter: String 'B' or None which decides whether you are looking (B)elow a certain threshold or not
        """
        
        data = data_setup()
        for i in range(len(data.ejec_data)):
            data.ex[i], data.ey[i], data.ez[i], data.vesc[i], data.ejec_KE[i], data.ejec_PE[i], data.Nclose[i], \
            data.ejected[i] = ejected_extract_final(data.IMBH_tracker[i], data.ejec_data[i], 'S') #Change S to E

        if filter == 'B':
            idx = np.where(data.ejec_KE < 1e46)                   # Change to 1e45 when more data
            ejected_KE1 = np.asarray([i for i in data.ejec_KE[idx]])
            ejected_PE1 = np.asarray([i for i in data.ejec_PE[idx]])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(data.ex, data.ey, data.ez)])[idx]

            idx = np.where(abs(ejected_PE1) < 1e46)
            ejected_KE = np.asarray([i/max(ejected_KE1[idx]) for i in ejected_KE1[idx]])
            ejected_PE = np.asarray([i/max(ejected_KE1[idx]) for i in ejected_PE1[idx]])
            final_pos = np.asarray(final_pos[idx])

        else:
            ejected_KE = np.asarray([i/max(data.ejec_KE) for i in data.ejec_KE])
            ejected_PE = np.asarray([i/max(data.ejec_KE) for i in data.ejec_PE])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(data.ex, data.ey, data.ez)])

        final_pos[np.isnan(final_pos)] = 0
        linex = np.linspace(-0.01*np.min(ejected_KE), 10)
        liney = -1*linex

        return ejected_KE, ejected_PE, final_pos, linex, liney

    def energy_plot_setup(self, KE, PE, linex, liney, cdata, clabel, no_axis, save_file, hist):
        """
        Function to plot the required data
        
        Inputs:
        KE/PE:       The x, y data of the plot
        linex/liney: The line of best fit to distinguish bound from unbound
        cdata:       The data for which the colour code is based off
        no_axis:     The number of axis plotting (zoom or not)
        save_file:   String to ensure no overwritten plots
        hist:        String to state whether plotting histogram or other coloured variable
        """

        plot_ini = plotter_setup()
        plt.figure(figsize=(10, 4))
        if no_axis == 2:
            gs = gridspec.GridSpec(1, 2)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0,1])
            ax1.set_ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
           
            for ax_ in [ax1, ax2]:
                ax_.plot(linex, liney, color = 'black', linestyle = '-.')
                ax_.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                plot_ini.tickers(ax_)
                
            if hist == 'Y':
                bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(KE, PE, bins=(200, 200), range=([0,1.1],[1.1*min(PE),0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE), 0]

                contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                #ax2.scatter(KE, PE, s = 15, color = 'red', edgecolors = 'black')
                #ax2.set_ylim(m(PE),0)
                ax2.set_xlim(0,max(KE))
                bin2d_sim, xedges_s, yedges_s, image = ax1.hist2d(KE, PE, bins=(80, 80), range=([0,0.35],[-0.35,0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 0.3, -0.3,0]
                contours = ax1.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                #ax1.scatter(KE, PE, s = 15, color = 'red', edgecolors = 'black')
                cbar = plt.colorbar(contours, ax=ax2)
                cbar.set_label(label = r'$\log_{10}(N/N_{max})$',rotation=270,labelpad=15)

                ax1.set_xlim(0,0.3)
                ax1.set_ylim(-0.3,0)

            if hist == 'N':
                colour_axes = ax1.scatter(KE, PE, c=cdata, s=2)
                colour_axes = ax2.scatter(KE, PE, c=cdata, s=2)
                plt.colorbar(colour_axes, ax=ax2, label = clabel)
                ax2.set_ylim(-1.05,0)
                ax2.set_xlim(0,1.05)
                ax1.set_xlim(0,0.3)
                ax1.set_ylim(-0.3,0)
                
            plt.savefig('figures/energy_diagram_'+str(save_file)+'_.pdf', dpi=500, bbox_inches='tight')
            return
        
        else:
            if hist == 'N':
                print(PE)
                fig, ax = plt.subplots()
                ax.set_title('Ejected Particles \nKinetic Energy vs. Potential Energy')
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                colour_axes = ax.scatter(KE, PE, c = cdata, s= 2)
                plt.colorbar(colour_axes, ax = ax, label = r'Final Distance to Core [pc]')
                plt.plot(linex, liney, color = 'black', linestyle = '-.')
                plt.xlim(1e-5,1.05)
                plt.ylim(-10000, -10**-15)
                #plt.yscale('symlog')
                #plt.xscale('log')
                plt.xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                plt.ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
                plot_ini.tickers(ax)
                plt.savefig('figures/energy_diagram_'+str(save_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
            else: 
                fig, ax = plt.subplots()

                bin2d_sim, xedges_s, yedges_s, image = ax.hist2d(KE, PE, bins=(100,300), range=([np.min(KE), 2],[1.1*min(PE),0.1]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE),0]

                ax.set_title('Ejected Particles \nKinetic Energy vs. Potential Energy')
                ax.set_xlim(1e-6, 1.1)
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                ax.set_ylim(-1, -10**-6)
                ax.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                ax.set_ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
                #ax.set_yscale('symlog')
                #ax.set_xscale('log')
                contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                #ax.scatter(KE, PE, s=15, color = 'red', edgecolors ='black')
                cbar = plt.colorbar(contours, ax=ax)
                cbar.set_label(r"$\log_{10}(N/N_{max})$", rotation=270,labelpad=15)
                
                ax.set_yscale('symlog')
                plot_ini.tickers(ax)
                plt.savefig('figures/energy_diagram_'+str(save_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
            return

    def KEPE_plotter(self):
            """
            Function to plot the KE vs. PE plot coloured with final positions or with histogram
            """

            save_file = ['fpos', 'fpos_crop']
            save_file_hist = ['histogram', 'histogram_crop']
            data_filt = [None, 'B']
            for i in range(2):
                axis = i+1
                ejected_KE, ejected_PE, final_pos, unbounded_x, unbounded_y = self.energy_data_extract(data_filt[i])
                self.energy_plot_setup(ejected_KE, ejected_PE, unbounded_x, unbounded_y, final_pos, 'Final Distance to Core [pc]', axis, save_file[i], 'N')  
                #self.finalpos_energy_plotter(ejected_KE, ejected_PE, unbounded_x, unbounded_y, final_pos, r'$\log_{10}n$', i, save_file_hist[i], 'Y')  

class vejection(object):
    """
    Class to plot the velocity ejection statistics
    """
    def plotter_func(self, xdata, ydata, cdata, clabel):
        """
        Function to plot the wanted data
        
        Inputs:
        xdata, ydata: The x and y variables wanting to plot
        cdata:        The data to be used as the colour code
        clabel:       Label describing the colour plot values
        """

        plot_ini = plotter_setup()
        fig, ax = plt.subplots()
        colour_axes = ax.scatter(xdata, ydata, c = cdata, zorder = 3)
        plt.colorbar(colour_axes, ax=ax, label = clabel)
        plot_ini.tickers(ax)
        ax.set_xlabel(r'Total IMBH Mass [$\frac{M}{10^3 M_{\odot}}$]')
        ax.set_ylabel(r'Ejection Velocity [km/s]')

        return ax

    def vejec_mean_plotter(self):
        """
        Function to plot how the total population/ mass of the system influences the ejection velocity
        """

        data = data_setup()
        for i in range(len(data.ejec_data)):
            data.ex[i], data.ey[i], data.ez[i], data.vesc[i], data.ejec_KE[i], data.ejec_PE[i], data.Nclose[i], \
            data.ejected[i] = ejected_extract_final(data.IMBH_tracker[i], data.ejec_data[i], 'S') #Change S to E
            
            data.vals_df = data.ejec_data[i].iloc[0]          
            data.tot_mass[i] = np.sum(data.vals_df[8].value_in(units.MSun))
            data.tot_pop[i] = data.vals_df[6]
            data.surv_time[i] = data.vals_df[-2].value_in(units.Myr)

        in_pop = np.unique(data.tot_pop)
        avg_vel = np.empty(len(in_pop))
        avg_surv = np.empty((len(in_pop)))

        iter = -1
        for pop_ in in_pop:
            iter += 1
            indices = np.where((data.tot_pop == pop_))[0]
            temp_escv = [data.vesc[i] for i in indices]
            temp_surv = [data.surv_time[i] for i in indices]

            avg_vel[iter] = np.mean(temp_escv)
            avg_surv[iter] = np.mean(temp_surv)

        avg_vel = np.asarray(avg_vel)
        #yint, slope = polyfit(popul, avg_vel, 1)
        #xlist = np.linspace(2.5,10.5)
        #ylist = [yint+slope*i for i in xlist]

        ax = self.plotter_func(in_pop, avg_vel, avg_surv, r'Average Ejection Time [Myr]')
        #ax.plot(xlist, ylist, linestyle = '--', zorder = 1)
        #ax.text(3, 50, r'$v_{ejec} = '+str('{:.2f}'.format(slope))+r'$N +'+str('{:.2f}'.format(yint)))
        ax.set_ylim(0, 600)
        ax.set_xlim(2.5, 10.5)
        plt.savefig('figures/mean_vej.pdf', dpi=300, bbox_inches='tight')

    def vejec_histogram(self):
        """
        Function to plot the ejection velocity w.r.t number of times occuring in bins
        """

        data = data_setup()

        for i in range(len(data.ejec_data)):
            data.ex[i], data.ey[i], data.ez[i], data.vesc[i], data.ejec_KE[i], data.ejec_PE[i], \
            data.Nclose[i], data.ejected[i] = ejected_extract_final(data.IMBH_tracker[i], data.ejec_data[i], 'S') #TO CHANGE WHEN MORE DATA
            
            data.vals_df = data.ejec_data[i].iloc[0]          
            data.tot_mass[i] = np.sum(data.vals_df[8].value_in(units.MSun))
            data.tot_pop[i] = data.vals_df[6]
            data.surv_time[i] = data.vals_df[-2].value_in(units.Myr)

        in_pop = np.unique(data.tot_pop)
        iter = -1
        colors = colour_picker()
        plt.title('Ejection Velocity Histogram')
        for pop_ in in_pop:
            if pop_ == 10 or pop_ == 50 or pop_ == 100:
                iter += 1
                indices = np.where((pop_ == in_pop))[0]
                esc_vel = data.vesc[indices[0]]
                n, bins, patches = plt.hist(esc_vel, 50, density=True, histtype = 'step', color=colors[iter])
            else:
                pass
        plt.xlabel(r'$v_{ejec}$ [km/s]')
        plt.ylabel(r'Fractional Occurence')
        plt.savefig('figures/vejection_histogram.pdf', dpi = 300, bbox_inches='tight')

class event_tracker(object):
    """
    Class to take stats of ejection vs. mergers
    """
    def __init__(self):
        data = data_setup()
        plot_init = plotter_setup()
                
        ejected = np.empty(len(data.ejec_data))
        for i in range(len(data.ejec_data)):
            data.ex[i], data.ey[i], data.ez[i], data.vesc[i], data.ejec_KE[i], data.ejec_PE[i], data.Nclose[i], \
            ejected[i] = ejected_extract_final(data.IMBH_tracker[i], data.ejec_data[i], 'A')

            data.vals_df = data.ejec_data[i].iloc[0] 
            data.tot_pop[i] = data.vals_df[6]

        in_pop = np.unique(data.tot_pop)
        frac_merge = np.empty(len(in_pop))
        
        iter = -1
        for pop_ in in_pop:
            iter +=1
            indices = np.where((data.tot_pop == pop_))[0]
            temp_frac = [ejected[i] for i in indices]
            frac_merge[iter] = np.mean(temp_frac)

        for i in range(len(data.ejec_data_GRX)):
            data.ex_GRX[i], data.ey_GRX[i], data.ez_GRX[i], data.vesc_GRX, data.ejec_KE_GRX[i], data.ejec_PE_GRX[i], \
            data.Nclose_GRX[i], data.ejected_GRX[i] = ejected_extract_final(data.IMBH_tracker_GRX[i], data.ejec_data_GRX[i], 'A')

            data.vals_df_GRX = data.ejec_data_GRX[i].iloc[0] 
            data.tot_pop_GRX[i] = data.vals_df_GRX[6]

        in_pop_GRX = np.unique(data.tot_pop_GRX)
        frac_merge_GRX = np.empty(len(in_pop_GRX))
        
        iter = -1
        for pop_ in in_pop_GRX:
            iter +=1
            indices = np.where((data.tot_pop_GRX == pop_))[0]
            temp_frac = [data.ejected_GRX[i] for i in indices]
            frac_merge_GRX[iter] = np.mean(temp_frac)


        fig, ax = plt.subplots()
        plot_init.tickers(ax)
        ax.set_xlabel(r'Number of IMBH [$N$]')
        ax.set_ylabel(r'Fraction of Simulations with SMBH Mergers')
        ax.set_ylim(0,1.05)
        ax.scatter(in_pop, frac_merge, color = 'red', edgecolors = 'black', label = 'Hermite')
        ax.scatter(in_pop_GRX, frac_merge_GRX, color = 'red', edgecolors = 'black', label = 'GRX')
        xints = [i for i in range(1+int(max(in_pop))) if i %10 == 0]
        ax.set_xticks(xints)
        plt.savefig('figures/SMBH_merge_fraction.pdf', dpi=300, bbox_inches='tight')

